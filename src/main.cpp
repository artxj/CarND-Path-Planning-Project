#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double METERS_AHEAD = 30.0;
const double LANE_WIDTH = 4.0;
const double TARGET_X = 30.0;
const unsigned int STEPS_COUNT = 50;
constexpr double TIME_STEP = 0.02;
const double CARS_MIN_DISTANCE = 30.0;
const double MAX_VELOCITY = 49.0;
const double VELOCITY_STEP = 0.25;

constexpr double time_per_velocity() { return TIME_STEP / 2.24; }

typedef struct {
  bool is_close;
  double velocity;
} CarInfo;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Reads sensor fusion information and provides an information about closest car in the same lane
void findClosestCar(const vector<vector<double> > &sensor_fusion,
  const double car_first_s, const double car_last_s, const int lane,
  const int old_lane, const unsigned int prev_size, bool &car_close, double &closest_car_velocity,
  double &closest_car_distance) {
  car_close = false;
  double min_distance = std::numeric_limits<double>::max();

  for (unsigned int i = 0; i < sensor_fusion.size(); ++i) {
    const vector<double> sensor_car = sensor_fusion[i];
    const double sensor_car_d = sensor_car[6];

    if (sensor_car_d > LANE_WIDTH * lane && sensor_car_d < LANE_WIDTH * (lane + 1)) {
      // car in same lane
      const double sensor_car_vx = sensor_car[3];
      const double sensor_car_vy = sensor_car[4];
      const double sensor_car_first_s = sensor_car[5];

      const double sensor_car_speed = sqrt(sensor_car_vx * sensor_car_vx +
        sensor_car_vy * sensor_car_vy);
      const double sensor_car_last_s = sensor_car_first_s + TIME_STEP * sensor_car_speed * prev_size;

      const double distance = sensor_car_last_s - car_last_s;
      if (distance > 0) {
        if (distance < CARS_MIN_DISTANCE) car_close = true;
        if (distance < min_distance) {
          // this car is closer than previous ones
          min_distance = distance;
          closest_car_velocity = sensor_car_speed;
        }
      }

      if ((sensor_car_first_s > car_first_s - 5) && (sensor_car_last_s < car_last_s + 5)) {
        car_close = true;
        min_distance = 0.0;
      }

    }
  }
  closest_car_distance = min_distance;
}

void calculatePath(const vector<double> &previous_path_x, const vector<double> &previous_path_y,
  const int lane, const double ref_velocity, const double car_s,
  const double car_x, const double car_y, const double car_yaw,
  const vector<double> &map_waypoints_s, const vector<double> &map_waypoints_x,
  const vector<double> &map_waypoints_y,
  vector<double> &next_x_vals, vector<double> &next_y_vals) {

  const unsigned int prev_size = previous_path_x.size();
  vector<double> ptsx;
  vector<double> ptsy;

  double ref_x = car_x;
  double ref_y = car_y;
  double ref_yaw = deg2rad(car_yaw);

  if (prev_size < 2) {
    const double prev_car_x = car_x - cos(ref_yaw);
    const double prev_car_y = car_y - sin(ref_yaw);

    ptsx.push_back(prev_car_x);
    ptsx.push_back(car_x);
    ptsy.push_back(prev_car_y);
    ptsy.push_back(car_y);
  } else {
    ref_x = previous_path_x[prev_size - 1];
    ref_y = previous_path_y[prev_size - 1];
    const double prev_ref_x = previous_path_x[prev_size - 2];
    const double prev_ref_y = previous_path_y[prev_size - 2];
    ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);

    ptsx.push_back(prev_ref_x);
    ptsx.push_back(ref_x);
    ptsy.push_back(prev_ref_y);
    ptsy.push_back(ref_y);
  }

  for (unsigned int i = 0; i < 3; ++i) {
    const vector<double> wp = getXY(car_s + METERS_AHEAD * (i + 1),
      LANE_WIDTH * (lane + 0.5),
      map_waypoints_s, map_waypoints_x, map_waypoints_y);
    ptsx.push_back(wp[0]);
    ptsy.push_back(wp[1]);
  }

  // shifting coords
  const double cos_ref_yaw = cos(-ref_yaw);
  const double sin_ref_yaw = sin(-ref_yaw);
  for (unsigned int i = 0; i < ptsx.size(); ++i) {
    const double shift_x = ptsx[i] - ref_x;
    const double shift_y = ptsy[i] - ref_y;
    ptsx[i] = shift_x * cos_ref_yaw - shift_y * sin_ref_yaw;
    ptsy[i] = shift_x * sin_ref_yaw + shift_y * cos_ref_yaw;
  }

  tk::spline sp;
  sp.set_points(ptsx, ptsy);

  for (unsigned int i = 0; i < prev_size; ++i) {
    next_x_vals.push_back(previous_path_x[i]);
    next_y_vals.push_back(previous_path_y[i]);
  }

  const vector<double> target = { TARGET_X, sp(TARGET_X) };
  const double target_dist = sqrt(target[0] * target[0] + target[1] * target[1]);

  double x_add = 0.0;
  const double n = target_dist / (time_per_velocity() * ref_velocity);

  for (unsigned int i = 0; i < STEPS_COUNT - prev_size; ++i) {
    double x = x_add + target[0] / n;
    double y = sp(x);
    x_add = x;

    double x_ref = x;
    double y_ref = y;

    x = x_ref * cos_ref_yaw + y_ref * sin_ref_yaw;
    y = x_ref * (-sin_ref_yaw) + y_ref * cos_ref_yaw;

    x += ref_x;
    y += ref_y;

    next_x_vals.push_back(x);
    next_y_vals.push_back(y);
  }
}

double calculateCost(const int lane, const double closest_car_distance,
  const double closest_car_velocity, const double ref_velocity) {
  if (closest_car_distance <= 0) return 1.0;
  double cost = 0.0;

  double velocity_diff = ref_velocity - closest_car_velocity;
  if (velocity_diff < 0) velocity_diff = 0;

  if (closest_car_distance < 3 * CARS_MIN_DISTANCE) {
    cost += 1 - exp(-velocity_diff / closest_car_distance);
  }

  if (lane != 1) cost += 0.2;

  if (cost > 1) cost = 1;
  return cost;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // new variables
  int lane = 1; // initial lane
  double ref_velocity = 0.0; // start with zero velocity

  h.onMessage([&lane, &ref_velocity,
    &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
    &map_waypoints_dx,&map_waypoints_dy]
    (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];
        	// Previous path's end s and d values
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

        	json msgJson;

          const unsigned int prev_size = previous_path_x.size();
          //cout << prev_size << endl;

          vector<double> costs = { 2, 2, 2 };
          vector<CarInfo> associated_info(3); // contains vectors of {is_car_close, velocity}
          for (unsigned int check_lane = 0; check_lane < 3; ++check_lane) {
            vector<double> x_vals;
            vector<double> y_vals;

            calculatePath(previous_path_x, previous_path_y, check_lane, ref_velocity,
              car_s, car_x, car_y, car_yaw, map_waypoints_s, map_waypoints_x,
              map_waypoints_y, x_vals, y_vals);

            const unsigned int next_size = x_vals.size();
            const double last_x = x_vals[next_size - 1];
            const double last_y = y_vals[next_size - 1];
            const double prev_last_x = x_vals[next_size - 2];
            const double prev_last_y = y_vals[next_size - 2];
            const double last_yaw = atan2(last_y - prev_last_y, last_x - prev_last_x);
            const vector<double> frenet = getFrenet(last_x, last_y, last_yaw,
              map_waypoints_x, map_waypoints_y);

            bool car_close = false;
            double closest_car_velocity = ref_velocity;
            double closest_car_distance = std::numeric_limits<double>::max();
            findClosestCar(sensor_fusion, car_s, frenet[0], check_lane, check_lane, prev_size,
              car_close, closest_car_velocity, closest_car_distance);

            costs[check_lane] = calculateCost(check_lane, closest_car_distance,
              closest_car_velocity, MAX_VELOCITY);
            //cout << "Cost for lane " << check_lane << " is " << costs[check_lane] << endl;
            CarInfo car_info;
            car_info.is_close = car_close;
            car_info.velocity = closest_car_velocity;

            associated_info[check_lane] = car_info;
          }
          if (lane == 0) costs[2] = 1;
          if (lane == 2) costs[0] = 1;

          const unsigned int prev_lane = lane;
          auto min_result = std::min_element(costs.begin(), costs.end());
          lane = std::distance(costs.begin(), min_result);

          const bool car_close = associated_info[lane].is_close;
          const double closest_car_velocity = associated_info[lane].velocity;
          if (car_close) {
            if (ref_velocity > closest_car_velocity) ref_velocity -= VELOCITY_STEP;
          } else if (ref_velocity < MAX_VELOCITY && prev_lane == lane) {
            ref_velocity += VELOCITY_STEP;
          }
          if (lane != prev_lane && ref_velocity > MAX_VELOCITY - 5 * VELOCITY_STEP && !car_close) {
            //ref_velocity -= 2 * VELOCITY_STEP;
          }

          /*
          if (prev_lane == 1 && associated_info[1].is_close) {
            lane = 0;
          }

          if (associated_info[1].is_close) {
            ref_velocity -= VELOCITY_STEP;
          } else if (ref_velocity < MAX_VELOCITY) {
            ref_velocity += VELOCITY_STEP;
          }
          */

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          calculatePath(previous_path_x, previous_path_y, lane, ref_velocity,
            car_s, car_x, car_y, car_yaw, map_waypoints_s, map_waypoints_x,
            map_waypoints_y, next_x_vals, next_y_vals);

        	// define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
      	  msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
