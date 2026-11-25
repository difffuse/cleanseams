// SPDX-License-Identifier: GPL-3.0-or-later
// Processes are spawned with an argument vector (posix_spawn), never through a shell.
#include "cleanseams/process.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/io.hpp"
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <poll.h>
#include <pthread.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
namespace cleanseams {
namespace {
[[noreturn]] void host(const std::string& field,const std::string& message,ErrorCode code=ErrorCode::host_failure) {throw ContractError(code,field,message);}
struct Pipe {
    int fd[2]{-1,-1};
    Pipe() {if(::pipe(fd)!=0)host("process","cannot create a pipe");}
    ~Pipe() {for(int& f:fd)if(f>=0)::close(f);}
    void close(int i) {if(fd[i]>=0) {::close(fd[i]);fd[i]=-1;}}
};
bool executable_file(const std::filesystem::path& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path,error) && ::access(path.c_str(),X_OK)==0;
}
}
ProcessResult run_process(const std::vector<std::string>& argv,std::string_view input,const ProcessOptions& options) {
    if(argv.empty())host("process","empty argument vector",ErrorCode::invalid_input);
    for(const auto& arg:argv)if(arg.find('\0')!=std::string::npos)host("process","embedded NUL argument",ErrorCode::invalid_input);
    ProcessResult result;
    if(options.cancel.requested()){result.cancelled=true;return result;}
    if(options.timeout.count()<0)host("process","negative timeout",ErrorCode::invalid_input);
    if(options.timeout.count()==0){result.timed_out=true;return result;}
    if(!executable_file(argv[0]))host("process","not an executable file: "+argv[0]);
    Pipe in,out,err;
    for(int f:{in.fd[0],in.fd[1],out.fd[0],out.fd[1],err.fd[0],err.fd[1]})
        if(::fcntl(f,F_SETFD,FD_CLOEXEC)<0)host("process","cannot isolate process pipes");
    for(int f:{in.fd[1],out.fd[0],err.fd[0]})
        if(::fcntl(f,F_SETFL,::fcntl(f,F_GETFL)|O_NONBLOCK)<0)host("process","cannot make process pipe nonblocking");
#ifdef F_SETNOSIGPIPE
    if(::fcntl(in.fd[1],F_SETNOSIGPIPE,1)<0)host("process","cannot suppress pipe signal");
#endif
    struct Spawn {
        posix_spawn_file_actions_t actions;posix_spawnattr_t attr;
        Spawn() {
            if(posix_spawn_file_actions_init(&actions))host("process","cannot initialize spawn actions");
            if(posix_spawnattr_init(&attr)){posix_spawn_file_actions_destroy(&actions);host("process","cannot initialize spawn attributes");}
        }
        ~Spawn(){posix_spawn_file_actions_destroy(&actions);posix_spawnattr_destroy(&attr);}
    } spawn;
    auto check=[](int code){if(code)host("process",std::string("cannot configure process: ")+std::strerror(code));};
    check(posix_spawn_file_actions_adddup2(&spawn.actions,in.fd[0],0));
    check(posix_spawn_file_actions_adddup2(&spawn.actions,out.fd[1],1));
    check(posix_spawn_file_actions_adddup2(&spawn.actions,err.fd[1],2));
    for(int f:{in.fd[0],in.fd[1],out.fd[0],out.fd[1],err.fd[0],err.fd[1]})check(posix_spawn_file_actions_addclose(&spawn.actions,f));
    check(posix_spawnattr_setpgroup(&spawn.attr,0));check(posix_spawnattr_setflags(&spawn.attr,POSIX_SPAWN_SETPGROUP));
    std::vector<char*> args;for(const auto& arg:argv)args.push_back(const_cast<char*>(arg.c_str()));args.push_back(nullptr);
    struct Child {
        pid_t pid=0;bool reaped=false;
        void stop(){if(pid>0)::kill(-pid,SIGKILL);}
        int reap(){int status=0;pid_t done;do{done=::waitpid(pid,&status,0);}while(done<0 && errno==EINTR);reaped=done==pid;return reaped?status:-1;}
        ~Child(){if(pid>0 && !reaped){stop();reap();}}
    } child;
    pid_t spawned_pid=0;
    check(posix_spawn(&spawned_pid,argv[0].c_str(),&spawn.actions,&spawn.attr,args.data(),environ));
    child.pid=spawned_pid;
    in.close(0);out.close(1);err.close(1);
#ifndef F_SETNOSIGPIPE
    // Block only this thread's SIGPIPE. Never change the embedding application's
    // process-wide handler while another concurrent runner may still be using it.
    struct PipeSignal {
        sigset_t old{},blocked{};bool pending=false;
        PipeSignal(){sigemptyset(&blocked);sigaddset(&blocked,SIGPIPE);const int e=pthread_sigmask(SIG_BLOCK,&blocked,&old);if(e)host("process","cannot block pipe signal");sigset_t p;sigpending(&p);pending=sigismember(&p,SIGPIPE)==1;}
        ~PipeSignal(){sigset_t p;sigpending(&p);if(!pending && sigismember(&p,SIGPIPE)==1){int signal;sigwait(&blocked,&signal);}pthread_sigmask(SIG_SETMASK,&old,nullptr);}
    } signal;
#endif
    const auto started=MonotonicClock::now();std::size_t written=0;
    if(input.empty())in.close(1);
    auto drain=[&](Pipe& pipe,std::string& text,bool& truncated) {
        char buffer[65536];
        // Bound each turn even if the child floods continuously, so cancellation
        // and deadlines cannot be starved by an always-readable pipe.
        for(int batch=0;batch<16 && pipe.fd[0]>=0;++batch) {
            const auto n=::read(pipe.fd[0],buffer,sizeof buffer);
            if(n>0){const auto size=static_cast<std::size_t>(n),keep=std::min(size,options.max_output-text.size());text.append(buffer,keep);truncated=truncated || keep<size;}
            else if(n==0){pipe.close(0);return;}
            else if(errno==EINTR)continue;
            else if(errno==EAGAIN || errno==EWOULDBLOCK)return;
            else host("process","cannot read process output");
        }
    };
    for(;;) {
        result.cancelled=options.cancel.requested();
        result.timed_out=!result.cancelled && MonotonicClock::now()-started>=options.timeout;
        // WNOWAIT keeps the direct PID reserved until its group has been stopped.
        siginfo_t info{};int status;
        do{status=::waitid(P_PID,static_cast<id_t>(child.pid),&info,WEXITED|WNOHANG|WNOWAIT);}while(status<0 && errno==EINTR);
        if(status<0)host("process","cannot inspect child status");
        if(info.si_pid || result.cancelled || result.timed_out) {
            child.stop();in.close(1);const int exit=child.reap();
            drain(out,result.out,result.out_truncated);drain(err,result.err,result.err_truncated);
            result.exit_code=exit<0?-1:WIFEXITED(exit)?WEXITSTATUS(exit):WIFSIGNALED(exit)?128+WTERMSIG(exit):-1;
            return result;
        }
        std::vector<pollfd> fds;
        if(in.fd[1]>=0)fds.push_back({in.fd[1],POLLOUT,0});
        if(out.fd[0]>=0)fds.push_back({out.fd[0],POLLIN,0});
        if(err.fd[0]>=0)fds.push_back({err.fd[0],POLLIN,0});
        const int ready=::poll(fds.data(),static_cast<nfds_t>(fds.size()),10);
        if(ready<0){if(errno==EINTR)continue;host("process","cannot poll process pipes");}
        for(const auto& p:fds) {
            if(!p.revents)continue;
            if(p.fd==in.fd[1]) {
                if(p.revents&(POLLERR|POLLHUP)){in.close(1);continue;}
                const auto n=::write(in.fd[1],input.data()+written,std::min<std::size_t>(input.size()-written,65536));
                if(n>0)written+=static_cast<std::size_t>(n);
                else if(n<0 && errno==EPIPE)in.close(1);
                else if(n<0 && errno!=EAGAIN && errno!=EINTR)host("process","cannot write process input");
                if(written==input.size())in.close(1);
            } else if(p.fd==out.fd[0])drain(out,result.out,result.out_truncated);
            else if(p.fd==err.fd[0])drain(err,result.err,result.err_truncated);
        }
    }
}
}
