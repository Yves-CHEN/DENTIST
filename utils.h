
#ifndef    __UTILS__
#define    __UTILS__
#endif   //__UTILS__


#include <string> 
using namespace std;
#include "headers.h"


void stop( const char *fmt, ...) 
{
    va_list args;
    va_start(args, fmt);
    vprintf (fmt, args);
    va_end(args);
    exit(-1);
}
//void stop(char* msg)
//{
//}

void to_upper(string &str)
{
    for(int i=0; i<str.size(); i++){
        if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
    }
}



void to_upper(char* str, int len)
{
    for(int i=0; i<len; i++){
        if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
    }
}


int split_string(const string &str, vector<string> &vec_str, string separator)
{
    if(str.empty()) return 0;
    vec_str.clear();

    int i=0;
    bool look=false;
    string str_buf;
    string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
    string::size_type pos;

    for(i=0; i<separator.size(); i++){
        pos=symbol_pool.find(separator[i]);
        if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
    }

    for(i=0; i<str.size(); i++){
        if( symbol_pool.find(str[i])!=string::npos )
        {
            if(!look) look=true;
            str_buf += str[i];
        }
        else
        {
            if(look){
                look=false;
                vec_str.push_back(str_buf);
                str_buf.erase(str_buf.begin(), str_buf.end());
            }
        }
    }
    if(look) vec_str.push_back(str_buf);

    return (int)vec_str.size();
}


