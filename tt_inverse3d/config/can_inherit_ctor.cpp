struct base {
    explicit base(int i) {}
};

struct derived : public base {
    using base::base;
};

int main() {
    derived d(42);
    return 0;
}
