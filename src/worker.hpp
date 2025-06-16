#pragma once

#include <memory>

#include <vector>
#include <atomic>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>

class Worker final
{
private:
    std::atomic<bool> active;
    std::atomic<bool> working;

    std::unique_ptr<std::thread> thread;
    std::mutex mutex;
    std::condition_variable condition;
    std::condition_variable condition_join;
    std::function<void()> job;

    void Loop();

public:
    Worker(std::function<void()>&& job);
    Worker(const Worker&) = delete;

    virtual ~Worker();

    void Terminate();
    void Notify();

    void Join();
};

class WorkerPool final
{
private:
    size_t max_workers_;
    std::vector<std::unique_ptr<Worker>> workers_;

    void Notify();
    void Join();

public:
    WorkerPool();
    virtual ~WorkerPool();

    void AddWorker(std::unique_ptr<Worker> worker);
    void Resolve();
    void Terminate();

    const size_t MaxWorkers() const
    {
        return max_workers_;
    }
};