#include <mphys/mphys.hpp>
#include <print>

int main()
{
    std::println("mphys v{}.{}", mphys::version_major, mphys::version_minor);
    return 0;
}
