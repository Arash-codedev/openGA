#pragma once

#include <chrono>

namespace OpenGA {
    template <typename ClockType>
    class StopWatch {
        typename ClockType::time_point start;

      public:
        StopWatch() : start(ClockType::now()) {}

        auto getDuration() {
            return std::chrono::duration_cast<std::chrono::milliseconds>(
                ClockType::now() - start);
        }
    };
} // namespace OpenGA
