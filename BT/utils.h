#include <iostream>
// ANSI color codes
#define RESET "\033[0m"
#define BLACK "\033[30m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"

// Custom color stream
class ColorStream
{
public:
    ColorStream(std::ostream &stream, const char *color)
        : m_stream(stream), m_color(color) {}

    template <typename T>
    ColorStream &operator<<(const T &val)
    {
        m_stream << m_color << val << RESET;
        return *this;
    }

    // Specialization for std::endl
    ColorStream &operator<<(std::ostream &(*pf)(std::ostream &))
    {
        pf(m_stream);
        return *this;
    }

private:
    std::ostream &m_stream;
    const char *m_color;
};

// Define color macros
#define black ColorStream(std::cout, BLACK)
#define red ColorStream(std::cout, RED)
#define green ColorStream(std::cout, GREEN)
#define yellow ColorStream(std::cout, YELLOW)
#define blue ColorStream(std::cout, BLUE)
#define magenta ColorStream(std::cout, MAGENTA)
#define cyan ColorStream(std::cout, CYAN)
#define white ColorStream(std::cout, WHITE)
