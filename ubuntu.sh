make distclean
sh configure --prefix=/usr
make -j8
sudo checkinstall --pkgname=real --pkgversion 0.0.31 --backup=no --default --deldoc
