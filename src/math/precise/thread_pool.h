// ============================================================================
// 线程池实现
// ============================================================================
//
// 提供高效的任务并行执行能力：
// - 预创建固定数量的工作线程
// - 任务队列支持动态添加任务
// - 支持同步等待所有任务完成
// - 线程安全的任务提交和执行
//
// 用于 NTT 并行计算等需要多线程的场景，避免频繁创建/销毁线程的开销。
// ============================================================================

#ifndef PRECISE_THREAD_POOL_H
#define PRECISE_THREAD_POOL_H

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace precise {

/**
 * @class ThreadPool
 * @brief 高效的线程池实现
 *
 * 使用方法：
 *   ThreadPool& pool = ThreadPool::instance();
 *   pool.submit([]() { do_work(); });
 *   pool.wait_all();  // 等待所有任务完成
 */
class ThreadPool {
public:
    /**
     * @brief 获取全局线程池实例（懒加载）
     * @return 线程池实例引用
     */
    static ThreadPool& instance() {
        static ThreadPool pool;
        return pool;
    }

    /**
     * @brief 提交任务到线程池
     * @param task 要执行的任务函数
     * @return future 对象，可用于等待任务完成
     */
    std::future<void> submit(std::function<void()> task) {
        std::packaged_task<void()> packaged_task(std::move(task));
        std::future<void> future = packaged_task.get_future();

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            tasks_.push(std::move(packaged_task));
        }

        condition_.notify_one();
        return future;
    }

    /**
     * @brief 提交多个任务并等待全部完成
     * @param tasks 任务列表
     *
     * 比逐个 submit + wait 更高效，减少同步开销。
     */
    void submit_and_wait(std::vector<std::function<void()>> tasks) {
        if (tasks.empty()) return;

        std::atomic<int> remaining{static_cast<int>(tasks.size())};
        std::mutex done_mutex;
        std::condition_variable done_condition;

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            for (auto& task : tasks) {
                std::packaged_task<void()> packaged_task([&, t = std::move(task)]() {
                    t();
                    if (--remaining == 0) {
                        std::lock_guard<std::mutex> lk(done_mutex);
                        done_condition.notify_one();
                    }
                });
                tasks_.push(std::move(packaged_task));
            }
        }

        condition_.notify_all();

        // 等待所有任务完成
        std::unique_lock<std::mutex> lock(done_mutex);
        done_condition.wait(lock, [&]() { return remaining == 0; });
    }

    /**
     * @brief 等待所有已提交任务完成
     *
     * 注意：此方法会阻塞直到任务队列为空且所有任务执行完毕。
     */
    void wait_all() {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        finished_condition_.wait(lock, [&]() {
            return tasks_.empty() && active_count_ == 0;
        });
    }

    /**
     * @brief 获取线程数量
     * @return 工作线程数量
     */
    size_t size() const { return workers_.size(); }

    /**
     * @brief 析构函数，停止所有线程
     */
    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            stop_ = true;
        }
        condition_.notify_all();

        for (auto& worker : workers_) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

private:
    /**
     * @brief 构造函数，创建工作线程
     *
     * 线程数量默认为硬件并发数 - 1（留一个给主线程）
     */
    ThreadPool() : stop_(false), active_count_(0) {
        size_t num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;
        if (num_threads > 1) num_threads--;  // 留一个给主线程
        if (num_threads < 2) num_threads = 2;  // 至少需要 2 个线程用于 NTT 并行

        for (size_t i = 0; i < num_threads; ++i) {
            workers_.emplace_back([this]() {
                while (true) {
                    std::packaged_task<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(queue_mutex_);
                        condition_.wait(lock, [&]() {
                            return stop_ || !tasks_.empty();
                        });

                        if (stop_ && tasks_.empty()) {
                            return;
                        }

                        task = std::move(tasks_.front());
                        tasks_.pop();
                        ++active_count_;
                    }

                    task();

                    {
                        std::lock_guard<std::mutex> lock(queue_mutex_);
                        --active_count_;
                        if (tasks_.empty() && active_count_ == 0) {
                            finished_condition_.notify_one();
                        }
                    }
                }
            });
        }
    }

    // 禁止拷贝
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    std::vector<std::thread> workers_;
    std::queue<std::packaged_task<void()>> tasks_;

    std::mutex queue_mutex_;
    std::condition_variable condition_;
    std::condition_variable finished_condition_;

    std::atomic<bool> stop_;
    std::atomic<int> active_count_;
};

} // namespace precise

#endif // PRECISE_THREAD_POOL_H