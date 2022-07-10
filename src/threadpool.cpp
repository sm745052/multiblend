#pragma once

#include "threadpool.h"
#include <thread>

Threadpool* Threadpool::instance;

/**********************************************************************
* Constructor (private)
**********************************************************************/
Threadpool::Threadpool(int _threads) {
	n_threads = _threads > 0 ? (std::min)((unsigned int)_threads, std::thread::hardware_concurrency()) : std::thread::hardware_concurrency();
	threads = new tp_struct[n_threads];

	for (int i = 0; i < n_threads; ++i) {
#ifdef _WIN32
		threads[i].handle = CreateThread(NULL, 1, (LPTHREAD_START_ROUTINE)Thread, &threads[i], 0, NULL);
#else
		pthread_create(&threads[i].handle, NULL, TP_Thread, &threads[i]);
#endif
		threads[i].main_mutex = &main_mutex;
		threads[i].return_mutex = &return_mutex;
		threads[i].main_cond = &main_cond;
		threads[i].return_cond = &return_cond;
		threads[i].free = true;
		threads[i].stop = false;
		threads[i].queue = &queue;
		threads[i].i = i;
	}
}

/**********************************************************************
* Destructor
**********************************************************************/
Threadpool::~Threadpool() {
	int i;

	{
		std::lock_guard<std::mutex> mlock(main_mutex);
		for (i = 0; i < n_threads; ++i) {
			threads[i].stop = true;
		}
	}

	main_cond.notify_all();
	for (i = 0; i < n_threads; ++i) {
#ifdef _WIN32
		WaitForSingleObject(threads[i].handle, INFINITE);
#else
		pthread_join(threads[i].handle, NULL);
#endif
	}
}

/**********************************************************************
* Threads
**********************************************************************/
#define P ((Threadpool::tp_struct*)param)

#ifdef _WIN32
DWORD WINAPI Threadpool::Thread(void* param) {
#else
void* TP_Thread(void* param) {
#endif
	while (true) {
		{
			std::unique_lock<std::mutex> mlock(*P->main_mutex);
			P->main_cond->wait(mlock, [=]{ return P->queue->size() || P->stop; });
			if (P->queue->size()) {
				P->function = P->queue->front();
				P->queue->pop_front();
				P->free = false;
			}
		}
		if (P->stop) break;

		P->function();

		{
			std::lock_guard<std::mutex> mlock(*P->return_mutex); // necessary
			P->free = true;
		}
		P->return_cond->notify_all();
	}

	return 0;
}

/**********************************************************************
* Wait
**********************************************************************/
void Threadpool::Wait() {
	if (!queue.size()) {
		int i;
		for (i = 0; i < n_threads; ++i) {
			if (!threads[i].free) break;
		}
		if (i == n_threads) return;
	}

	{
		std::unique_lock<std::mutex> rlock(return_mutex);
		return_cond.wait(rlock, [=]{
			if (queue.size()) return false;
			for (int i = 0; i < n_threads; ++i) {
				if (!threads[i].free) {
					return false;
				}
			}
			return true;
		});
	}
}

/**********************************************************************
* Queue
**********************************************************************/
void Threadpool::Queue(std::function<void()> function) {
	std::lock_guard<std::mutex> mlock(main_mutex); // not sure what this is for
	queue.push_back(std::move(function));
	main_cond.notify_one(); // changed from notify_all()
}
