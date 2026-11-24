

#ifndef _PROGRESS_BAR_
#define _PROGRESS_BAR_

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <atomic>
#include <thread>

class ProgressBar {
  public:
  ProgressBar(int total, std::ostream& os = std::cout, int width = 60)
      : total_(total), width_(width), os_(os),
        start_(std::chrono::steady_clock::now()) {}

  void update(int done) {
    using namespace std::chrono;
    double elapsed = duration<double>(steady_clock::now() - start_).count();
    double frac    = (total_ > 0) ? (double)done / total_ : 0.0;
    double eta     = (frac > 0.05) ? elapsed / frac * (1.0 - frac) : -1.0;

    int filled = (int)(frac * width_);
    std::string bar(filled, '=');
    if (filled < width_) { bar += '>'; bar += std::string(width_ - filled - 1, ' '); }

    if (done < total_){
      os_ << '\r'
        << '[' << std::left << std::setw(width_) << bar << ']'
        << ' ' << std::right << std::setw(3) << (int)(frac * 100) << '%'
        // << "  elapsed " << fmt_time(elapsed)
        << "  eta " << (eta >= 0 ? fmt_time(eta) : "?")
        // << "  (" << done << '/' << total_ << ')'
        << "   "
        << std::flush;
    }else{
      os_ << '\r'
        << '[' << std::left << std::setw(width_) << bar << ']'
        << ' ' << std::right << std::setw(3) << (int)(frac * 100) << '%'
        << "  total " << fmt_time(elapsed)
        // << "  eta " << (eta >= 0 ? fmt_time(eta) : "?")
        // << "  (" << done << '/' << total_ << ')'
        << endl;
    }    
  }

  private:
  int total_, width_;
  std::ostream& os_;
  std::chrono::steady_clock::time_point start_;

  static std::string fmt_time(double s) {
      std::ostringstream ss;
      int h   = (int)s / 3600;
      int m   = ((int)s % 3600) / 60;
      int sec = (int)s % 60;
      if (h > 0)
          ss << h << 'h' << std::setfill('0') << std::setw(2) << m
             << 'm' << std::setw(2) << sec << 's';
      else if (m > 0)
          ss << m << 'm' << std::setfill('0') << std::setw(2) << sec << 's';
      else
          ss << std::fixed << std::setprecision(1) << s << 's';
      return ss.str();
  }
};




#endif