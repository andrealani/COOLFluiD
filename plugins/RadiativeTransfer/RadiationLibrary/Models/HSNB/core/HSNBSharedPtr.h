#ifndef COOLFluiD_RadiativeTransfer_SHARED_PTR_H
#define COOLFluiD_RadiativeTransfer_SHARED_PTR_H

// Use std::shared_ptr implementation available in C++11 if supported
#if __cplusplus > 199711L
    //#pragma message("using std::shared_ptr")
    #include <memory>
    template <typename T> using HSNBSharedPtr = std::shared_ptr<T>;

// Otherwise, use our own implementation (probably slower)
#else

#include <algorithm> // for std::swap


/**
 * Very simple reference counting smart pointer.
 */
template <typename T>
class HSNBSharedPtr
{
public:
    HSNBSharedPtr() : mp_ref(NULL) { }

    /**
     * Constructs a new HSNBSharedPtr owning the object pointed to by the pointer.
     */
    explicit HSNBSharedPtr(T* ptr) : mp_ref(new RefCount(ptr)) { }

    /**
     * Copy constructor.
     */
    HSNBSharedPtr(const HSNBSharedPtr<T>& other) :
        mp_ref(other.mp_ref)
    {
        if (mp_ref != NULL)
            mp_ref->count++;
    }

    /**
     * Decreases reference count to the object.  If this is the last reference,
     * then the object is deleted.
     */
    ~HSNBSharedPtr()
    {
        // Protect against non-initialized pointers
        if (mp_ref == NULL) return;

        // Reduce reference count and delete object
        if (--(mp_ref->count) == 0) {
            delete mp_ref->ptr;
            delete mp_ref;
        }
    }

    /// Assignment operator.
    HSNBSharedPtr<T>& operator = (HSNBSharedPtr<T> other) {
        std::swap(mp_ref, other.mp_ref);
        return *this;
    }

    /// Equality operator.
    bool operator == (HSNBSharedPtr<T> other) {
        return (mp_ref == other.mp_ref);
    }

    /// Returns reference to object owned by this shared pointer.
    T& operator*() const { return *(mp_ref->ptr); };

    /// Returns pointer to object owned by this shared pointer.
    T* operator->() const { return mp_ref->ptr; };

    template <class U, class V>
    friend bool operator == (const HSNBSharedPtr<U>&, const HSNBSharedPtr<V>&);

    template <class U, class V>
    friend bool operator != (const HSNBSharedPtr<U>&, const HSNBSharedPtr<V>&);

private:

    // Helper class for managing the reference count.
    struct RefCount
    {
        RefCount() : ptr(NULL), count(0) { }
        RefCount(T* p) : ptr(p), count(1) { }
        T*  ptr;
        int count;
    };

private:

    RefCount* mp_ref;

}; // class HSNBSharedPtr

// Compare two HSNBSharedPtr objects.
template <class T, class U>
bool operator == (const HSNBSharedPtr<T>& lhs, const HSNBSharedPtr<U>& rhs) {
    return (lhs.mp_ref == rhs.mp_ref);
}
template <class T, class U>
bool operator != (const HSNBSharedPtr<T>& lhs, const HSNBSharedPtr<U>& rhs) {
    return (lhs.mp_ref != rhs.mp_ref);
}

#endif // __cplusplus > 199711L
#endif // SHARED_PTR_H

