CPPFLAGS := -g -O2 -Wall -Wno-unused-function -Wextra -Wno-missing-field-initializers
jFloaty.exe: jFloaty.cpp *.h
	g++ ${CPPFLAGS} $< -o $@

clean:
	rm -f jFloaty.exe
