#include <iostream>
#include <time.h>
#include <unistd.h>

using namespace std;

// Timer class for simple measurement of execution times
class Timer {
	private:
		clock_t start, finish;
	
	public:
		void tic() {
			start=time(NULL);
		}

		clock_t toc() {
			finish = time(NULL);

			return (  (finish - start));
		}

};
