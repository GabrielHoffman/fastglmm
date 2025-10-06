

#ifndef _PROGRESS_BAR_
#define _PROGRESS_BAR_

#include <iostream>
#include <chrono>
#include <thread> // For std::this_thread::sleep_for


class ProgressBar {

  public:

  ProgressBar( const int &total ):
    total(total) {
    startTime = std::chrono::high_resolution_clock::now();
  }

  // Update with number of completed tasks and out-stream
  void update( const int &progress, ostream &strm);

  private:
  std::chrono::high_resolution_clock::time_point startTime;
  int total;
  int barWidth = 50; // Width of the progress bar in characters
};


void ProgressBar::update( const int &progress, ostream &strm){

  // Calculate progress percentage
  double percentage = static_cast<double>(progress) / total;
  int filledWidth = static_cast<int>(barWidth * percentage);

  // Print the bar
  strm << "\r["; // \r returns cursor to beginning of line
  for (int i = 0; i < filledWidth; ++i) {
      strm << "=";
  }
  for (int i = filledWidth; i < barWidth; ++i) {
      strm << " ";
  }
  strm << "] " << static_cast<int>(percentage * 100.0) << "%";

  if( progress > 0 ){
    // Calculate ETA
    auto currentTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

    int minutes, seconds;

    if( progress == total){
      minutes = static_cast<int>(elapsedTime / 60);
      seconds = static_cast<int>(elapsedTime) % 60;
    }else{
      double estimatedTotalTime = (elapsedTime / percentage);
      double etaSeconds = estimatedTotalTime - elapsedTime;

      minutes = static_cast<int>(etaSeconds / 60);
      seconds = static_cast<int>(etaSeconds) % 60;
    }
    strm << " ETA: " << minutes << "m " << seconds << "s";
    strm << std::flush; // Ensure immediate output
  }
}




#endif