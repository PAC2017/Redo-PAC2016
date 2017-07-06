协方差算法优化

1. 程序介绍：
程序提供4个源文件，main.f90, check.f90, Covariance.f90和Covariance1.f90，以及一个Makefile文件。
其中：

a) main.f90提供数据初始化，调用功能计算函数，结果验证。

b) check.f90辅助校验。且check.f90文件IsMissingPheno函数中的判断条件不会改变。
c) Covariance1.f90与Covariance.f90均是功能计算函数，目前函数逻辑实现一致，用于验证使用。


直接make可以生成可执行文件 main，直接运行可执行文件得出程序计时。
编译：
$ make
ifort -c main.f90
ifort -c Covariance1.f90
ifort -c Covariance.f90
ifort -c check.f90
ifort  main.o Covariance1.o Covariance.o check.o -o main

运行：
./main
 total time(s):    运行时间    
 verification: correct

2. 比赛规则：
a) 不得修改 main.f90, check.f90 和 Covariance1.f90 三个文件，若修改，本环节成绩直接作废。
b) 本测试程序在调用Covariance函数时，参数传递值只是是一种常见情况，
选手需要根据程序逻辑进行优化，不得只针对传递的参数值做特定优化，否则本环节成绩作废。
c) 程序结果验证通过，成绩有效，否则本环节成绩作废。
d) 比赛提供测试平台KNL，建议使用平台预安装的Intel compiler 2016.0.0.3.
e) 最终优化比较，以程序输出计时为评判依据，用时越短者成绩越高。


