* 仅提供xming-fonts下载
    * 搜索后得到xming下载地址

* 在WSL2中连接被拒绝
    * 命令应改为export DISPLAY=$(route.exe print | grep 0.0.0.0 | head -1 | awk '{print $4}'):0.0
    * https://theunixtips.com/xming-client-4-rejected-from-ip/ 将IP加入Xming的hosts配置文件

* mesh的绘制出现问题，构造函数调用错误，因为没有理解程序结构，通过查看mesh.cpp得到方法