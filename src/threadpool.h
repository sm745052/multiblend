#pragma once

#include <mutex>
#include <condition_variable>
#include <functional>
#include <deque>
#ifdef _WIN32
#include <Windows.h>
#endif

#ifndef _WIN32
void* TP_Thread(void* param);
#endif

class Threadpool {
public:
	static Threadpool* GetInstance(int threads = 0) { if (!instance) instance = new Threadpool(threads); return instance; }
	void Queue(std::function<void()> function);
	int GetNThreads() { return n_threads; };
	void Wait();
	struct tp_struct {
#ifdef _WIN32
		HANDLE handle;
#else
		pthread_t handle;
#endif
		std::function<void()> function;
		bool free;
		bool stop;
		std::mutex* main_mutex;
		std::mutex* return_mutex;
		std::condition_variable* main_cond;
		std::condition_variable* return_cond;
		std::deque<std::function<void()>>* queue;
		int i;
	};

private:
	static Threadpool* instance;
	Threadpool(int _threads = 0); // constructor is private
	~Threadpool();
#ifdef _WIN32
	static DWORD WINAPI Thread(void* param);
#endif
	tp_struct* threads;
	std::deque<std::function<void()>> queue;
	int n_threads;
	std::mutex main_mutex;
	std::mutex return_mutex;
	std::condition_variable main_cond;
	std::condition_variable return_cond;
};
