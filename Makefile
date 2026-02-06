CPPFLAGS := -O2
jFloaty.exe: jFloaty.cpp *.h
	g++ ${CPPFLAGS} $< -o $@
