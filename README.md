# libtsp

An exact & heuristic TSP solver written in C++.

## Example usage

The following example showcases the libtsp api:
```cpp
    TspInputGraphEdge edges[2] = {
        {.cost = 10.0, .from = 1, .to = 2},
        {.cost = 10.0, .from = 2, .to = 1}
    };
    constexpr TspInputGraphDescriptor graph{
        .edges = edges,
        .num_edges = 2
    };
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    tspAsymmetricSolve(&graph, &optionsDesc, &outputDesc);
```

## Prerequisites

- Git
- CMake (3.22.1 or higher)
- C++ compiler with C++20 support (MSVC 17+, gcc 11+ or clang 12+)

## Supported architectures

libtsp aims to be compatible with all architectures.
Feel free to create issues for architecture-induced compilation failures.

## Building

### Installing prerequisites

In this section we propose a method of installing the required prerequisites for building libtsp on Windows, macOS and Linux.

#### Windows

With the winget package manager installed & up-to-date from the Microsoft Store, you can install the prerequisites as
follows:

```bash
winget install Microsoft.VisualStudio.2022.Community --silent --override "--wait --quiet --add ProductLang En-us --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"
winget install Git.Git --source winget
winget install Kitware.CMake --source winget
```

After installing these packages, make sure to refresh your PATH by restarting your explorer.exe in the Task Manager and
opening a new Terminal launched by said explorer.

#### macOS

```bash
xcode-select --install # if not already installed

# install Homebrew package manager if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install git # if not already installed by xcode command line tools
brew install cmake
```

#### Ubuntu

```bash
sudo apt update
sudo apt install -y build-essential
sudo apt install -y git
sudo apt install -y cmake
```

### Building the native library & other targets

To build all native targets, run the following commands valid for both Windows with PowerShell and Unix-like systems
starting from the root directory of the repository:

```bash
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release --parallel
```

## Testing

### C++ Tests

To run the C++ unit tests, starting from the root directory of the repository, run the following commands valid for both
Windows with PowerShell and Unix-like systems:

```bash
cd build
ctest --verbose --build-config Release --output-on-failure
```


## License

This project is licensed under the MIT License.

## Contributing

Contributions are welcome! Please submit issues and pull requests to help improve libtsp.
