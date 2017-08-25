#include <sys/mman.h>

unsigned char *my_mmap(size_t len, int fd) {
        void *result = mmap(0, len, PROT_READ, MAP_SHARED, fd, 0);
        return (unsigned char*)( result == MAP_FAILED ? 0 : result );
}

void my_munmap(void *len, unsigned char *p) {
        munmap( p, (size_t)len ) ; 
}
