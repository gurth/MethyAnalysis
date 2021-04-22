## Build save plug-in

CMAKE_SOURCE_DIR_:=../../../../..

all: build

build: libsave.so
	cp libsave.so $(CMAKE_SOURCE_DIR_)/plug-in/

libsave.so: save.c ../../profile_node.h $(CMAKE_SOURCE_DIR_)/include/config.h
	g++ -std=c++11 -fPIC -shared -I $(CMAKE_SOURCE_DIR_)/include -o libsave.so save.c

clean:
	rm -rf libsave.so $(CMAKE_SOURCE_DIR_)/plug-in/libsave.so

help:
	@echo ""
	@echo "Usage:"
	@echo "\tmake target"
	@echo ""
	@echo "Targets:"
	@echo "\tall:\t\tThe default if no target is provided."
	@echo "\tbuild:\t\tBuild all."
	@echo "\tlibsave.so:\t\tBuild libsave.so."
	@echo "\thelp:\t\tShow this message."