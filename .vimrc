set nocompatible
if has('mouse')
  set mouse=a
endif

set autoindent
set shiftwidth=2
set tabstop=2
set expandtab

call plug#begin('~/.vim/plugged')

Plug 'lervag/vimtex'

call plug#end()
  
let g:vimtex_view_general_viewer = 'SumatraPDF.exe'
let g:vimtex_view_general_options = '-reuse-instance -forward-search @tex @line @pdf'
let g:vimtex_view_general_options_latexmk = '-reuse-instance'

let fortran_do_enddo=1
set showmatch
set cursorline
