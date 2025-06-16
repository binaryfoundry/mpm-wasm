#include "worker.hpp"

Worker::Worker(std::function<void()>&& job) :
    job(std::move(job)),
    active(true),
    working(false)
{
    thread = std::make_unique<std::thread>([this] {
        Loop();
    });
}

Worker::~Worker()
{
    Terminate();
}

void Worker::Terminate()
{
    if (active)
    {
        active = false;
        condition.notify_one();
        thread->join();
    }
}

void Worker::Notify()
{
    working = true;
    condition.notify_one();
}

void Worker::Join()
{
    if (working)
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition_join.wait(lock, [this] {
            return !working || !active;
        });
        working = false;
    }
}

void Worker::Loop()
{
    do
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this] {
            return working || !active;
        });

        if (active && working)
        {
            job();
        }

        condition_join.notify_one();
        working = false;
    } while (active);
}

WorkerPool::WorkerPool()
{
    max_workers_ = std::thread::hardware_concurrency();
}

WorkerPool::~WorkerPool()
{
    Terminate();
}

void WorkerPool::AddWorker(std::unique_ptr<Worker> worker)
{
    if (workers_.size() >= max_workers_)
        throw std::runtime_error("Max number of workers exceeded.");

    workers_.push_back(std::move(worker));
}

void WorkerPool::Resolve()
{
    Notify();
    Join();
}

void WorkerPool::Notify()
{
    for (auto& w : workers_)
    {
        w->Notify();
    }
}

void WorkerPool::Join()
{
    for (auto& w : workers_)
    {
        w->Join();
    }
}

void WorkerPool::Terminate()
{
    for (auto& w : workers_)
    {
        w->Terminate();
    }
}
