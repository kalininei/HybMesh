#include "core_hmconnection_java.h"
#include "stdlib.h"
#include "client.h"

void exception(JNIEnv* env, const char* err){
	jclass Exception = (*env)->FindClass(env, "java/lang/RuntimeException");
	(*env)->ThrowNew(env, Exception, err);
}

JNIEXPORT jint JNICALL Java_HybMesh_require_1connection(JNIEnv* env, jobject obj){
	const char* server_path;
	int err;
	jclass self;
	jfieldID fid;
	jstring pth;
	int ret = -1;

	/* get path to server */
	self = (*env)->GetObjectClass(env, obj);
	fid = (*env)->GetFieldID(env, self, "hmserverpath", "Ljava/lang/String;");
	if (fid == NULL) exception(env, "failed to find object field");
	pth = (*env)->GetObjectField(env, obj, fid);
	server_path = (*env)->GetStringUTFChars(env, pth, NULL);
	if (server_path == NULL) exception(env, "failed to get server application path");
	printf("called with %s\n", server_path);

	/* create connection */
	err = HybmeshClientToServer_new(server_path, &ret);
	(*env)->ReleaseStringUTFChars(env, pth, server_path);

	if (err < 0) exception(env, "failed to establish a connection");
	return (jint)ret;
}

JNIEXPORT jbyte JNICALL Java_HybMesh_get_1signal(JNIEnv* env, jobject obj, jint fd){
	int err;
	char sig = '0';
	err = HybmeshClientToServer_get_signal(fd, &sig);
	if (err < 0) exception(env, "failed to get a signal from the pipe");
	return (jbyte)sig;
}

JNIEXPORT jbyteArray JNICALL Java_HybMesh_get_1data(JNIEnv* env, jobject obj, jint fd){
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

JNIEXPORT void JNICALL Java_HybMesh_send_1signal(JNIEnv* env, jobject obj, jint fd, jbyte sig){
	int err = HybmeshClientToServer_send_signal(fd, sig);
	if (err < 0) exception(env, "failed to send a signal to the pipe");
}

JNIEXPORT void JNICALL Java_HybMesh_send_1data(
		JNIEnv* env, jobject obj, jint fd, jbyteArray data){
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

JNIEXPORT void JNICALL Java_HybMesh_break_1connection(JNIEnv* env, jobject obj, jint fd){
	int err = HybmeshClientToServer_delete(fd);
	if (err < 0) exception(env, "failed to break hybmesh connection");
}
