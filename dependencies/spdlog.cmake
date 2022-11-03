FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.9.2
)
set(SPDLOG_BUILD_SHARED ON)
FetchContent_MakeAvailable(spdlog)
