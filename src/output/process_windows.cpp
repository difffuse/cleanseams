// SPDX-License-Identifier: GPL-3.0-or-later
// Native Windows process ownership and lifetime management.
#include "cleanseams/process.hpp"
#include "cleanseams/errors.hpp"
#include <windows.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <limits>
#include <thread>
namespace cleanseams {
namespace {
[[noreturn]] void host(const std::string& message,ErrorCode code=ErrorCode::host_failure) {throw ContractError(code,"process",message);}
struct Handle {
    HANDLE value=nullptr;
    explicit Handle(HANDLE h=nullptr):value(h){}
    Handle(const Handle&)=delete;Handle& operator=(const Handle&)=delete;
    ~Handle(){close();}
    void close(){if(value && value!=INVALID_HANDLE_VALUE)CloseHandle(value);value=nullptr;}
};
std::wstring wide(const std::string& text) {
    if(text.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))host("argument too long",ErrorCode::invalid_input);
    if(text.empty())return {};
    const auto length=static_cast<int>(text.size());
    const int n=MultiByteToWideChar(CP_UTF8,MB_ERR_INVALID_CHARS,text.data(),length,nullptr,0);
    if(!n)host("invalid UTF-8 argument",ErrorCode::invalid_input);
    std::wstring out(static_cast<std::size_t>(n),L'\0');
    if(!MultiByteToWideChar(CP_UTF8,MB_ERR_INVALID_CHARS,text.data(),length,out.data(),n))host("cannot convert argument");
    return out;
}
std::wstring quoted(const std::string& text) {
    std::wstring out=L"\"";std::size_t slashes=0;
    for(wchar_t c:wide(text)) {
        if(c==L'\\'){++slashes;continue;}
        out.append(slashes*(c==L'"'?2:1),L'\\');slashes=0;
        if(c==L'"')out+=L'\\';out+=c;
    }
    out.append(slashes*2,L'\\');return out+L'"';
}
}
ProcessResult run_process(const std::vector<std::string>& argv,std::string_view input,const ProcessOptions& options) {
    if(argv.empty())host("empty argument vector",ErrorCode::invalid_input);
    for(const auto& arg:argv)if(arg.find('\0')!=std::string::npos)host("embedded NUL argument",ErrorCode::invalid_input);
    ProcessResult result;
    if(options.cancel.requested()){result.cancelled=true;return result;}
    if(options.timeout.count()<0)host("negative timeout",ErrorCode::invalid_input);
    if(options.timeout.count()==0){result.timed_out=true;return result;}
    std::wstring command;for(const auto& arg:argv){if(!command.empty())command+=L' ';command+=quoted(arg);}
    if(command.size()>=32767)host("command line exceeds Windows limit",ErrorCode::invalid_input);
    const auto application=wide(argv.front());
    Handle in_read,in_write,out_read,out_write,err_read,err_write;
    SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES),nullptr,TRUE};
    if(!CreatePipe(&in_read.value,&in_write.value,&security,0) || !CreatePipe(&out_read.value,&out_write.value,&security,0) ||
       !CreatePipe(&err_read.value,&err_write.value,&security,0))host("cannot create process pipes");
    for(HANDLE h:{in_write.value,out_read.value,err_read.value})if(!SetHandleInformation(h,HANDLE_FLAG_INHERIT,0))host("cannot isolate process pipes");
    Handle job(CreateJobObjectW(nullptr,nullptr));if(!job.value)host("cannot create process job");
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};limits.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    if(!SetInformationJobObject(job.value,JobObjectExtendedLimitInformation,&limits,sizeof(limits)))host("cannot set process job limits");
    // A handle allowlist prevents concurrent sessions/Blender handles leaking into
    // the child. The job is assigned atomically at creation (Windows 10+).
    SIZE_T bytes=0;InitializeProcThreadAttributeList(nullptr,2,0,&bytes);
    std::vector<unsigned char> storage(bytes);
    auto* attributes=reinterpret_cast<LPPROC_THREAD_ATTRIBUTE_LIST>(storage.data());
    if(!InitializeProcThreadAttributeList(attributes,2,0,&bytes))host("cannot initialize process attributes");
    struct AttributeGuard {LPPROC_THREAD_ATTRIBUTE_LIST p;~AttributeGuard(){DeleteProcThreadAttributeList(p);}} guard{attributes};
    HANDLE inherit[]{in_read.value,out_write.value,err_write.value};
    if(!UpdateProcThreadAttribute(attributes,0,PROC_THREAD_ATTRIBUTE_HANDLE_LIST,inherit,sizeof(inherit),nullptr,nullptr) ||
       !UpdateProcThreadAttribute(attributes,0,PROC_THREAD_ATTRIBUTE_JOB_LIST,&job.value,sizeof(job.value),nullptr,nullptr))host("cannot set process attributes");
    STARTUPINFOEXW startup{};startup.StartupInfo.cb=sizeof(startup);startup.StartupInfo.dwFlags=STARTF_USESTDHANDLES;
    startup.StartupInfo.hStdInput=in_read.value;startup.StartupInfo.hStdOutput=out_write.value;startup.StartupInfo.hStdError=err_write.value;
    startup.lpAttributeList=attributes;PROCESS_INFORMATION info{};
    if(!CreateProcessW(application.c_str(),command.data(),nullptr,nullptr,TRUE,CREATE_NO_WINDOW|EXTENDED_STARTUPINFO_PRESENT,
                       nullptr,nullptr,&startup.StartupInfo,&info))host("cannot start executable (Windows error "+std::to_string(GetLastError())+")");
    Handle process(info.hProcess),main_thread(info.hThread);in_read.close();out_write.close();err_write.close();
    std::array<std::thread,3> workers;std::array<std::exception_ptr,3> failures{};std::atomic<bool> failed=false;
    auto protect=[&](std::size_t i,const auto& work){try{work();}catch(...){failures[i]=std::current_exception();failed.store(true);}};
    auto read=[&](HANDLE pipe,std::string& target,bool& truncated) {
        char buffer[8192];DWORD count=0;
        for(;;) {
            if(!ReadFile(pipe,buffer,sizeof(buffer),&count,nullptr)) {
                const auto error=GetLastError();if(error==ERROR_BROKEN_PIPE || error==ERROR_OPERATION_ABORTED)return;
                host("cannot read process output (Windows error "+std::to_string(error)+")");
            }
            if(!count)return;
            const auto keep=std::min<std::size_t>(count,options.max_output-target.size());
            target.append(buffer,keep);truncated=truncated || keep<count;
        }
    };
    auto stop=[&]{job.close();WaitForSingleObject(process.value,INFINITE);};
    auto join=[&]{for(auto& worker:workers)if(worker.joinable())worker.join();};
    try {
        workers[0]=std::thread([&]{protect(0,[&]{read(out_read.value,result.out,result.out_truncated);});});
        workers[1]=std::thread([&]{protect(1,[&]{read(err_read.value,result.err,result.err_truncated);});});
        workers[2]=std::thread([&]{protect(2,[&]{
            std::size_t offset=0;
            while(offset<input.size()) {
                DWORD sent=0;const auto count=static_cast<DWORD>(std::min<std::size_t>(input.size()-offset,8192));
                if(!WriteFile(in_write.value,input.data()+offset,count,&sent,nullptr)) {
                    const auto error=GetLastError();if(error==ERROR_BROKEN_PIPE || error==ERROR_NO_DATA || error==ERROR_OPERATION_ABORTED)break;
                    host("cannot write process input (Windows error "+std::to_string(error)+")");
                }
                if(!sent)break;offset+=sent;
            }
        });in_write.close();});
        const auto started=MonotonicClock::now();
        for(;;) {
            result.cancelled=options.cancel.requested();
            result.timed_out=!result.cancelled && MonotonicClock::now()-started>=options.timeout;
            const auto status=WaitForSingleObject(process.value,0);
            if(status==WAIT_FAILED)host("cannot wait for executable");
            if(status==WAIT_OBJECT_0 || result.cancelled || result.timed_out || failed.load())break;
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        // Kill ordinary descendants even after a successful direct-child exit;
        // they must not keep stdin/stdout handles alive and hang the joins.
        stop();join();
        for(const auto& error:failures)if(error)std::rethrow_exception(error);
        DWORD exit=0;if(!GetExitCodeProcess(process.value,&exit))host("cannot read exit status");
        result.exit_code=static_cast<int>(exit);return result;
    } catch(...){stop();join();throw;}
}
}
