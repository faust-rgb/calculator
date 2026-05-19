// ============================================================================
// BigIntData 实现 (SSO - Small Object Optimization)
// ============================================================================

#include "math/precise/precise_decimal.h"
#include <algorithm>

namespace precise {

// ============================================================================
// BigIntData 构造函数和析构函数
// ============================================================================

BigIntData::BigIntData() : ptr(sso), size_(0), capacity_(SSO_CAP) { }

BigIntData::BigIntData(size_t n, uint32_t val) : ptr(sso), size_(static_cast<uint32_t>(n)), capacity_(SSO_CAP) {
    if (n > SSO_CAP) {
        ptr = new uint32_t[n];
        capacity_ = static_cast<uint32_t>(n);
    }
    std::fill(ptr, ptr + size_, val);
}

BigIntData::BigIntData(const BigIntData& other) : ptr(sso), size_(other.size_), capacity_(SSO_CAP) {
    if (other.size_ > SSO_CAP) {
        ptr = new uint32_t[other.size_];
        capacity_ = other.size_;
    }
    std::copy(other.ptr, other.ptr + other.size_, ptr);
}

BigIntData::BigIntData(BigIntData&& other) noexcept : ptr(sso), size_(other.size_), capacity_(SSO_CAP) {
    if (other.ptr != other.sso) {
        ptr = other.ptr;
        capacity_ = other.capacity_;
        other.ptr = other.sso;
        other.capacity_ = SSO_CAP;
    } else {
        std::copy(other.sso, other.sso + other.size_, sso);
    }
    other.size_ = 0;
}

BigIntData::BigIntData(std::initializer_list<uint32_t> list) : ptr(sso), size_(static_cast<uint32_t>(list.size())), capacity_(SSO_CAP) {
    if (list.size() > SSO_CAP) {
        ptr = new uint32_t[list.size()];
        capacity_ = static_cast<uint32_t>(list.size());
    }
    std::copy(list.begin(), list.end(), ptr);
}

BigIntData::BigIntData(const uint32_t* first, const uint32_t* last) : ptr(sso), size_(0), capacity_(SSO_CAP) {
    size_t n = last - first;
    if (n > SSO_CAP) {
        ptr = new uint32_t[n];
        capacity_ = static_cast<uint32_t>(n);
    }
    size_ = static_cast<uint32_t>(n);
    std::copy(first, last, ptr);
}

BigIntData::~BigIntData() {
    if (ptr != sso) delete[] ptr;
}

// ============================================================================
// BigIntData 赋值操作符
// ============================================================================

BigIntData& BigIntData::operator=(const BigIntData& other) {
    if (this == &other) return *this;
    if (other.size_ > capacity_) {
        if (ptr != sso) delete[] ptr;
        ptr = new uint32_t[other.size_];
        capacity_ = other.size_;
    }
    size_ = other.size_;
    std::copy(other.ptr, other.ptr + other.size_, ptr);
    return *this;
}

BigIntData& BigIntData::operator=(BigIntData&& other) noexcept {
    if (this == &other) return *this;
    if (ptr != sso) delete[] ptr;
    if (other.ptr != other.sso) {
        ptr = other.ptr;
        capacity_ = other.capacity_;
        other.ptr = other.sso;
        other.capacity_ = SSO_CAP;
    } else {
        ptr = sso;
        capacity_ = SSO_CAP;
        std::copy(other.sso, other.sso + other.size_, sso);
    }
    size_ = other.size_;
    other.size_ = 0;
    return *this;
}

// ============================================================================
// BigIntData 容量操作
// ============================================================================

void BigIntData::push_back(uint32_t val) {
    if (size_ == capacity_) reserve(capacity_ * 2);
    ptr[size_++] = val;
}

void BigIntData::reserve(size_t n) {
    if (n <= capacity_) return;
    uint32_t* new_ptr = new uint32_t[n];
    std::copy(ptr, ptr + size_, new_ptr);
    if (ptr != sso) delete[] ptr;
    ptr = new_ptr;
    capacity_ = static_cast<uint32_t>(n);
}

void BigIntData::resize(size_t n, uint32_t val) {
    if (n > capacity_) reserve(n);
    if (n > size_) std::fill(ptr + size_, ptr + n, val);
    size_ = static_cast<uint32_t>(n);
}

// ============================================================================
// BigIntData 插入和删除操作
// ============================================================================

void BigIntData::erase(const uint32_t* first, const uint32_t* last) {
    size_t count = last - first;
    size_t start_idx = first - ptr;
    if (count == 0) return;
    std::move(ptr + start_idx + count, ptr + size_, ptr + start_idx);
    size_ -= static_cast<uint32_t>(count);
}

void BigIntData::insert(uint32_t* pos, const uint32_t* first, const uint32_t* last) {
    size_t count = last - first;
    size_t offset = pos - ptr;
    if (size_ + count > capacity_) reserve(std::max(size_ + count, (size_t)capacity_ * 2));
    std::move_backward(ptr + offset, ptr + size_, ptr + size_ + count);
    std::copy(first, last, ptr + offset);
    size_ += static_cast<uint32_t>(count);
}

void BigIntData::insert(uint32_t* pos, size_t count, uint32_t val) {
    size_t offset = pos - ptr;
    if (size_ + count > capacity_) reserve(std::max(size_ + count, (size_t)capacity_ * 2));
    std::move_backward(ptr + offset, ptr + size_, ptr + size_ + count);
    std::fill(ptr + offset, ptr + offset + count, val);
    size_ += static_cast<uint32_t>(count);
}

} // namespace precise