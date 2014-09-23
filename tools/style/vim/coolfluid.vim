" place this in your $VIMHOME (.vim on unix)
" and source it from your vim config file (or by hand by using source)

filetype plugin indent on

autocmd FileType cpp set textwidth=78 shiftwidth=3 formatoptions+=tcroq 
autocmd FileType cpp set expandtab 



