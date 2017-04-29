# CarND UKF Project

My implementation for CarND UKF project (Term2-P2).

## Params

```cpp
std_a_ = 0.9;

std_yawdd_ = 0.55;

P_ = MatrixXd::Identity(n_x_, n_x_);
```

## Compile
Under `Ubuntu 16.04 LTS` with `gcc/g++ >= v5.4`.

```
mkdir build && cd build
cmake .. && make -j3
```

## Results

The output results all satisfy the rubrics.

```
./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../output.txt
```
```
RMSE
0.0649021
0.0828896
 0.330927
 0.212981
```
```
./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt ../data/output1.txt 
```
```
RMSE
0.0771659
0.0847319
 0.598939
 0.583213
```
```
./UnscentedKF ../data/sample-laser-radar-measurement-data-2.txt ../data/output2.txt 
```
```
RMSE
0.192247
0.190768
0.412487
0.491183
```