// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
namespace cleanseams::output_detail {
// Python float(text) over ASCII: optional surrounding whitespace (\t\n\v\f\r, \x1c-\x1f,
// space), an optional sign, then digits with an optional fraction and exponent, where
// single underscores may separate digits (PEP 515: "1_000.5", not "1_.5" or "1__0"), or
// inf/infinity/nan (case-insensitive). Hexadecimal and anything else is not a number.
inline std::optional<double> python_float_text(std::string_view text) {
    constexpr std::string_view space=" \t\n\r\f\v\x1c\x1d\x1e\x1f";
    const auto a=text.find_first_not_of(space);if(a==std::string_view::npos)return std::nullopt;
    const auto b=text.find_last_not_of(space);std::string body(text.substr(a,b-a+1));
    std::string lower=body;for(auto& ch:lower)ch=static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    std::size_t i=0;double sign=1;
    if(i<lower.size() && (lower[i]=='+' || lower[i]=='-')) {sign=lower[i]=='-'?-1:1;++i;}
    const auto rest=std::string_view(lower).substr(i);
    if(rest=="inf" || rest=="infinity")return sign*std::numeric_limits<double>::infinity();
    if(rest=="nan")return std::numeric_limits<double>::quiet_NaN();
    // digitpart: digit (["_"] digit)*; the digits are copied without underscores.
    std::string clean=lower.substr(0,i);
    auto digitpart=[&](std::size_t& j) {
        std::size_t n=0;
        while(j<lower.size()) {
            if(std::isdigit(static_cast<unsigned char>(lower[j]))) {clean.push_back(lower[j]);++j;++n;}
            else if(lower[j]=='_' && n && j+1<lower.size() && std::isdigit(static_cast<unsigned char>(lower[j+1])))++j;
            else break;
        }
        return n;
    };
    std::size_t j=i;std::size_t digits=digitpart(j);
    if(j<lower.size() && lower[j]=='.') {clean.push_back('.');++j;digits+=digitpart(j);}
    if(!digits)return std::nullopt;
    if(j<lower.size() && lower[j]=='e') {
        clean.push_back('e');++j;if(j<lower.size() && (lower[j]=='+' || lower[j]=='-'))clean.push_back(lower[j++]);
        if(!digitpart(j))return std::nullopt;
    }
    if(j!=lower.size())return std::nullopt;
    return std::strtod(clean.c_str(),nullptr);
}
}
