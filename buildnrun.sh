#!bin/bash

cmake --preset debug -Wno-dev

cp build/debug/compile_commands.json build/compile_commands.json

cmake --build --preset debug --target krylov

bin/debug/krylov
