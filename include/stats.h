#ifndef GRAPHPI_STATS_H
#define GRAPHPI_STATS_H

#include <string>
#include <vector>

#define STATS_TABLE \
X(MOVE_IDX, "move_idx") \
X(CMP_COND, "cmp_cond")   \
X(PUSH_VAL, "push_val")

#define TIMER_TABLE \
Y(TOTAL, "total") \
Y(INTERSECT, "intersect")

#define X(a, b) a,
enum StatsName {
    STATS_TABLE
};
#undef X

#define X(a, b) b,
constexpr const char *stats_names[] = {
        STATS_TABLE
};
constexpr unsigned STATS_NAME_SIZE = sizeof(stats_names) / sizeof(char *);
#undef X

#define Y(a, b) a,
enum TimerName {
    TIMER_TABLE
};
#undef Y

#define Y(a, b) b,
constexpr const char *timer_names[] = {
        TIMER_TABLE
};
constexpr unsigned TIMER_NAME_SIZE = sizeof(timer_names) / sizeof(char *);
#undef Y

constexpr bool COUNT = false;
constexpr bool TIMER = true;

struct Stats {
    Stats() : counters(STATS_NAME_SIZE),  timers(TIMER_NAME_SIZE) {}

    std::vector<long long unsigned> counters;

    std::vector<double> timers;

    void inc(StatsName name) {
        if constexpr (COUNT) {
            ++counters[name];
        }
    }

    void add(StatsName name, unsigned val) {
        if constexpr (COUNT) {
            counters[name] += val;
        }
    }

    void start(double& val) const ;

    void end(TimerName name, double start);

    std::string ToString() const {
        std::string str = "Stats{ ";
        if constexpr (COUNT) {
            for (unsigned i = 0; i < STATS_NAME_SIZE; ++i) {
                str += std::string(stats_names[i]) + ":" + std::to_string(counters[i]) + ", ";
            }
            str += "\b\b";
        }
        if constexpr (TIMER) {
            for(unsigned i = 0; i < TIMER_NAME_SIZE; ++i) {
                str += std::string(timer_names[i]) + ":" + std::to_string(timers[i]) + "s, ";
            }
            str += "\b\b";
        }
        str += " }";
        return str;
    }
};

#endif //GRAPHPI_STATS_H
