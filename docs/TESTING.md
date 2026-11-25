# Testing

The default build runs two native smoke checks with data constructed in memory.

```sh
cmake --preset native-debug
cmake --build --preset native-debug
ctest --preset native-debug
```

`api_mesh_analysis` calls the C++ library on a quad and checks its area,
triangulation, corner count and advertised features. `cli_roundtrip` checks
the CLI serialization of an empty mesh. Neither check needs a host application,
Python interpreter or external geometry.

The [SDK example](../examples/installed/CMakeLists.txt) compiles the same C++
consumer against an installed package. Follow the [build instructions](BUILD.md)
to configure it in a separate project.
