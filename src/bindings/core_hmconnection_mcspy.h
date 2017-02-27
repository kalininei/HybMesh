#ifndef CORE_HYBMESH_CONNECTION_CS_H
#define CORE_HYBMESH_CONNECTION_CS_H

#ifdef __cplusplus
extern "C" {
#endif

int require_connection(const char* path);
void send_signal(int id, char sig);
void send_data(int id, int sz, char* data);
char get_signal(int id);
int get_data(int id, char** data);
int get_data1(int id);
void get_data2(int id, int sz, char* data);
void break_connection(int id);
void free_char_array(char* data);


#ifdef __cplusplus
}
#endif

#endif
