FetchContent_Declare(
        cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11
        GIT_TAG v2.3.1
)
FetchContent_MakeAvailable(cli11)

target_include_directories(tool PUBLIC ${cli11_SOURCE_DIR}/include/)
