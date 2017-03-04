#include "core_hmconnection_java.h"
#include "stdlib.h"
#include "client.h"

void exception(JNIEnv* env, const char* err){
	jclass Exception = (*env)->FindClass(env, "java/lang/RuntimeException");
	(*env)->ThrowNew(env, Exception, err);
}

JNIEXPORT jint JNICALL Java_Hybmesh_00024Worker_require_1connection(
		JNIEnv *env, jobject obj, jstring str){
	const char* server_path;
	int err;
	jclass self;
	jfieldID fid;
	int ret = -1;

	/* get path to server */
	server_path = (*env)->GetStringUTFChars(env, str, NULL);
	if (server_path == NULL) exception(env, "failed to get server application path");

	/* create connection */
	err = HybmeshClientToServer_new(server_path, &ret);
	(*env)->ReleaseStringUTFChars(env, str, server_path);

	if (err < 0) exception(env, "failed to establish a connection");
	return (jint)ret;
}

JNIEXPORT jbyte JNICALL Java_Hybmesh_00024Worker_get_1signal(
		JNIEnv * env, jobject obj, jint fd){
	int err;
	char sig = '0';
	err = HybmeshClientToServer_get_signal(fd, &sig);
	if (err < 0) exception(env, "failed to get a signal from the pipe");
	return (jbyte)sig;
}

JNIEXPORT jbyteArray JNICALL Java_Hybmesh_00024Worker_get_1data(
		JNIEnv *env, jobject obj, jint fd){
	int err;
	int sz;
	char* data;
	jbyteArray ret;

	err = HybmeshClientToServer_get_data(fd, &sz, &data);
	if (err < 0) exception(env, "failed to get data from the pipe");

	/* copy to jbyteArray */
	ret = (*env)->NewByteArray(env, sz);
	if (ret == NULL) exception(env, "failed to allocate the buffer");
	(*env)->SetByteArrayRegion(env, ret, 0 , sz, data);

	/* free cside and return */
	free(data);
	return ret;
}

JNIEXPORT void JNICALL Java_Hybmesh_00024Worker_send_1signal(
		JNIEnv *env, jobject obj, jint fd, jbyte sig){
	int err = HybmeshClientToServer_send_signal(fd, sig);
	if (err < 0) exception(env, "failed to send a signal to the pipe");
}

JNIEXPORT void JNICALL Java_Hybmesh_00024Worker_send_1data(
		JNIEnv * env, jobject obj, jint fd, jbyteArray data){
	char* cdata;
	int sz;
	int err;

	cdata = (*env)->GetByteArrayElements(env, data, NULL);
	if (cdata == NULL) exception(env, "failed to get buffer");
	sz = (*env)->GetArrayLength(env, data);

	err = HybmeshClientToServer_send_data(fd, sz, cdata);
	if (err < 0) exception(env, "failed to send data to the pipe");

	(*env)->ReleaseByteArrayElements(env, data, cdata, 0);
}

JNIEXPORT void JNICALL Java_Hybmesh_00024Worker_break_1connection(
		JNIEnv *env, jobject obj, jint fd){
	int err = HybmeshClientToServer_delete(fd);
	if (err < 0) exception(env, "failed to break hybmesh connection");
}
