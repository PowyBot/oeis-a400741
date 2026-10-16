CFLAGS = -O3 -march=native -pthread -Wall -Wextra

all: orbit_count line_count

%: src/%.c src/gf2.h
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f orbit_count line_count
