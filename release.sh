#!/bin/bash

date_suffix=`date +%Y%m%d%H%M%S`
hash=`git rev-parse --short HEAD`

cd ..
SRC_PATH=echo-cancellation
RELEASE_PATH=AEC-lib
# 删除旧文件
rm -rf $RELEASE_PATH

# 创建文件夹并拷贝源文件
mkdir -p $RELEASE_PATH
cp -r $SRC_PATH $RELEASE_PATH/jni
cd $RELEASE_PATH/jni

# 编译
ndk-build

# 拷贝静态库
cd ..
cp -r obj/local/arm64-v8a/ .
cp -r obj/local/armeabi-v7a/ .

# 删除多余文件
rm -rf jni/.vscode/
rm -rf jni/.git/
rm -rf arm64-v8a/objs/
rm -rf armeabi-v7a/objs/
rm -rf obj/

# 打包
cd ..
tar zcvf $RELEASE_PATH-${date_suffix}-${hash}.tar.gz $RELEASE_PATH
rm -rf $RELEASE_PATH
