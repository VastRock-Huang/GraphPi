#include "stats.h"
#include "common.h"

void Stats::end(TimerName name, double start) {
    if constexpr (TIMER) {
        timers[name] += (get_wall_time() - start);
    }
}

void Stats::start(double &val) const {
    if constexpr (TIMER) {
        val = get_wall_time();
    }
}

