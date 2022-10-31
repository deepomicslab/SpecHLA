# To install Novoalign, you need a permission, so please go to the http://www.novocraft.com/products/novoalign/ ask for permission. And remember link the binary file of Novoalign to the bin/ folder. 
# if you meet error "libarpack.so.2: cannot open shared object file: No such file or directory" while using SpecHap, please try:
LD_LIBRARY_PATH=/usr/local/lib
sudo apt update
sudo apt install libarpack++2-dev
