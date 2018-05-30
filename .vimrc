"Enable line number 
set number

" Plugin configuration
call plug#begin()
Plug 'scrooloose/nerdtree' "Arbo
Plug 'Xuyuanp/nerdtree-git-plugin' "Git color arbo
Plug 'iamcco/markdown-preview.vim' "Markdown preview for vim
Plug 'iamcco/mathjax-support-for-mkdp' "Markfown math
Plug 'biosugar0/singularity-vim' "Singularity synthax hl
Plug 'lervag/vimtex' "Vim editor
" A Vim Plugin for Lively Previewing LaTeX PDF Output
Plug 'xuhdev/vim-latex-live-preview', { 'for': 'tex' }
" 3 plugs for vim themes and code recognition
Plug 'itchyny/lightline.vim'
Plug 'sheerun/vim-polyglot'
Plug 'joshdick/onedark.vim'
Plug 'broadinstitute/vim-wdl' "Wdl syntax highlight
"  ...
call plug#end()

filetype plugin indent on

"Nerdtree Conf
autocmd vimenter * NERDTree
autocmd StdinReadPre * let s:std_in=1
autocmd VimEnter * if argc() == 0 && !exists("s:std_in") | NERDTree | endif
autocmd StdinReadPre * let s:std_in=1
autocmd VimEnter * if argc() == 1 && isdirectory(argv()[0]) && !exists("s:std_in") | exe 'NERDTree' argv()[0] | wincmd p | ene | endif

"Auto-ident
set expandtab
set shiftwidth=2
set softtabstop=2
set autoindent

"Close parenthesis
inoremap {      {}<Left>
inoremap {<CR>  {<CR>}<Esc>O
inoremap {{     {
inoremap {}     {}

 " Split keyboard shortcuts
 nnoremap <silent> ¬ <c-w>l
 nnoremap <silent> Ì <c-w>h
 nnoremap <silent> È <c-w>k
 nnoremap <silent> Ï <c-w>j

"Link to pdf viewer
let g:livepreview_previewer = 'open -a Preview'
autocmd Filetype tex setl updatetime=1
nmap <F12> :LLPStartPreview<cr>

" Theme activation 
let g:onedark_termcolors=256
syntax on
colorscheme onedark
