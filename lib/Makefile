CC = gcc
CFLAGS = -Wall -shared -fPIC -O3
LDFLAGS = -lm

TARGET = segplot.so
SRC = segplot.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(TARGET)