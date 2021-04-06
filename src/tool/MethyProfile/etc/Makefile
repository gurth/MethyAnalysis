## Build save plug-in

all: build

build: libsave.so

libsave.so: save.c ../profile_node.h
	g++ -fPIC -shared -o libsave.so save.c

clean:
	rm -rf libsave.so

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