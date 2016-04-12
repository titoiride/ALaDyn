
extern "C" { 

#if defined (__xlC__)
 void memaddr(void * var,unsigned long long int * addr ) {
#elif defined (_MSC_VER)
 void MEMADDR(void * var,unsigned long long int * addr ) {
#else
 void memaddr_(void * var,unsigned long long int * addr ) {
#endif
    *addr =  (unsigned long long int)var;
    return;
  }

}


