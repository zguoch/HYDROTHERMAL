# 编译
直接在src目录下执行makefile，因为module有依赖关系，所以需要按照先后顺序执行
```bash
make modules_mach
make mod
make fmod
make obj
make program
```

