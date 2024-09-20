#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <unistd.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include "fmt/format.h"
#include "memusage.h"
#include "macros.h"

namespace {
#if defined(__linux__) || defined(__linux)
std::map<std::string, std::string> proc_fields(pid_t pid, std::string const & endpoint)
{
    std::map<std::string, std::string> D;
    const std::string str = fmt::format("/proc/{}/{}", pid, endpoint);
    std::ifstream f(str);
    for(std::string key ; f >> key ; ) {
        std::string value;
        ASSERT_ALWAYS(!key.empty() && key.back() == ':');
        key.erase(key.end() - 1);
        std::getline(f, value);
        D[key] = value;
    }
    return D;
}
#endif
#if defined(__APPLE__)
bool get_mach_task_info(task_t task, mach_task_basic_info_data_t * info)
{
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    kern_return_t ret = task_info(task, MACH_TASK_BASIC_INFO, 
		(task_info_t) info, &count);
    return ret == KERN_SUCCESS && count == MACH_TASK_BASIC_INFO_COUNT;
}
#endif
}


/* Returns memory usage, in KiB (1024 bytes).
 * This is the VmSize field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure. When not supported, simply return 0.
 */
size_t
Memusage (void)
{
#if defined(__linux__) || defined(__linux)
    size_t mem;
    if (std::istringstream(proc_fields(getpid(), "status")["VmSize"]) >> mem)
        return mem;
    else
        return -1;
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    if (get_mach_task_info(mach_task_self(), &info))
        /* it seems that virtual_size is often inaccessible */
        return info.resident_size >> 10U;
    else
        return -1;
#else
    return 0;
#endif
}

/* same as above, for resident memory (column RES of top) */
size_t
Memusage2 (void)
{
#if defined(__linux__) || defined(__linux)
    size_t mem;
    if (std::istringstream(proc_fields(getpid(), "status")["VmRSS"]) >> mem)
        return mem;
    else
        return -1;
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    if (get_mach_task_info(mach_task_self(), &info))
        return info.resident_size >> 10U;
    else
        return -1;
#else
    return 0;
#endif
}

/* Returns peak memory usage, in KiB (1024 bytes).
 * This is the VmPeak field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
size_t
PeakMemusage (void)
{
#if defined(__linux__) || defined(__linux)
    size_t mem;
    if (std::istringstream(proc_fields(getpid(), "status")["VmPeak"]) >> mem)
        return mem;
    else
        return -1;
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    if (get_mach_task_info(mach_task_self(), &info))
        return info.resident_size_max >> 10U;
    else
        return -1;
#else
    return 0;
#endif
}
