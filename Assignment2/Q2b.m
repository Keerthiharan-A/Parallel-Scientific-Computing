clc
clear
serial_jacobi_1_x = [0 0.00437633 0.00865365 0.0128329 0.0169152 0.0209013 0.0247924 0.0285894 0.0322934 0.0359054 0.0394263 0.0428572 0.0461993 0.0494535 0.0526209 0.0557027 0.0586999 0.0616137 0.0644452 0.0671955 0.069866 0.0724577 0.074972 0.07741 0.0797729 0.0820622 0.084279 0.0864248 0.0885008 0.0905083 0.0924489 0.0943238 0.0961344 0.0978822 0.0995686 0.101195 0.102763 0.104274 0.105729 0.107131 0.108479 0.109777 0.111025 0.112225 0.113378 0.114486 0.115551 0.116573 0.117555 0.118497 0.119402 0.12027 0.121103 0.121903 0.12267 0.123407 0.124114 0.124793 0.125445 0.126071 0.126673 0.127252 0.127808 0.128343 0.128858 0.129354 0.129833 0.130294 0.130739 0.131168 0.131583 0.131985 0.132373 0.132749 0.133113 0.133466 0.133808 0.13414 0.134462 0.134775 0.135079 0.135374 0.13566 0.135938 0.136208 0.13647 0.136724 0.136969 0.137206 0.137436 0.137656 0.137869 0.138072 0.138267 0.138453 0.138629 0.138795 0.138952 0.139097 0.139232 0.139356 0.139467 0.139566 0.139653 0.139726 0.139786 0.139831 0.139861 0.139875 0.139874 0.139856 0.139821 0.139768 0.139697 0.139607 0.139498 0.139368 0.139218 0.139047 0.138854 0.138638 0.1384 0.138139 0.137853 0.137543 0.137208 0.136848 0.136462 0.136049 0.13561 0.135143 0.134648 0.134125 0.133573 0.132993 0.132382 0.131742 0.131071 0.13037 0.129637 0.128873 0.128077 0.127248 0.126387 0.125492 0.124564 0.123601 0.122604 0.121572 0.120505 0.119402 0.118262 0.117086 0.115872 0.11462 0.11333 0.112 0.110631 0.109222 0.107772 0.106279 0.104745 0.103167 0.101546 0.0998793 0.0981672 0.0964086 0.0946025 0.092748 0.090844 0.0888894 0.0868833 0.0848244 0.0827117 0.0805438 0.0783196 0.0760377 0.0736968 0.0712956 0.0688327 0.0663066 0.0637157 0.0610587 0.058334 0.0555399 0.0526748 0.0497372 0.0467252 0.0436373 0.0404716 0.0372264 0.0339 0.0304904 0.0269959 0.0234147 0.0197448 0.0159844 0.0121316 0.00818457 0.00414133 0];
serial_jacobi_1_y = [0 0.0037684 0.00743857 0.0110123 0.0144913 0.0178774 0.0211723 0.0243779 0.0274958 0.0305278 0.0334757 0.0363412 0.0391261 0.041832 0.0444608 0.0470141 0.0494936 0.0519011 0.0542383 0.0565068 0.0587083 0.0608445 0.0629171 0.0649277 0.0668779 0.0687693 0.0706037 0.0723825 0.0741073 0.0757799 0.0774016 0.078974 0.0804988 0.0819773 0.0834111 0.0848018 0.0861507 0.0874593 0.0887291 0.0899615 0.0911578 0.0923196 0.0934482 0.0945448 0.095611 0.0966479 0.0976569 0.0986393 0.0995964 0.100529 0.101439 0.102328 0.103196 0.104045 0.104875 0.105689 0.106487 0.107269 0.108039 0.108795 0.109541 0.110275 0.111 0.111716 0.112425 0.113127 0.113824 0.114515 0.115203 0.115887 0.116569 0.11725 0.11793 0.11861 0.119292 0.119974 0.12066 0.121348 0.12204 0.122737 0.123438 0.124146 0.124859 0.12558 0.126308 0.127044 0.127789 0.128542 0.129305 0.130078 0.130862 0.131656 0.132462 0.133279 0.134109 0.134951 0.135805 0.136672 0.137553 0.138447 0.139356 0.140278 0.141214 0.142165 0.14313 0.144111 0.145106 0.146116 0.147141 0.148181 0.149237 0.150308 0.151394 0.152495 0.153612 0.154743 0.15589 0.157052 0.158229 0.15942 0.160627 0.161847 0.163083 0.164332 0.165595 0.166871 0.168161 0.169465 0.17078 0.172109 0.173449 0.174801 0.176165 0.177539 0.178924 0.180319 0.181724 0.183138 0.18456 0.18599 0.187428 0.188873 0.190324 0.191781 0.193243 0.19471 0.19618 0.197653 0.199129 0.200606 0.202084 0.203562 0.205039 0.206514 0.207988 0.209457 0.210923 0.212383 0.213837 0.215284 0.216724 0.218154 0.219573 0.220982 0.222379 0.223762 0.22513 0.226483 0.227819 0.229138 0.230437 0.231715 0.232972 0.234205 0.235415 0.236598 0.237755 0.238883 0.239981 0.241048 0.242082 0.243082 0.244046 0.244973 0.245861 0.246709 0.247514 0.248276 0.248992 0.249662 0.250283 0.250853 0.251371 0.251835 0.252243 0.252593 0.252884 0.253113 0.25328 0.25338 0.253414];
parallel_jacobi_2_y = [0 0.00380064 0.00750302 0.0111089 0.0146199 0.0180379 0.0213646 0.0246017 0.0277509 0.030814 0.0337927 0.0366887 0.0395038 0.0422396 0.0448978 0.0474801 0.0499883 0.0524239 0.0547888 0.0570845 0.0593127 0.0614751 0.0635733 0.0656089 0.0675836 0.0694989 0.0713565 0.0731579 0.0749048 0.0765987 0.0782411 0.0798337 0.0813778 0.0828751 0.084327 0.085735 0.0871006 0.0884252 0.0897104 0.0909574 0.0921678 0.0933429 0.0944841 0.0955927 0.0966702 0.0977178 0.0987368 0.0997286 0.100694 0.101635 0.102553 0.103448 0.104323 0.105177 0.106013 0.106831 0.107633 0.108419 0.109191 0.10995 0.110697 0.111433 0.112159 0.112876 0.113585 0.114286 0.114982 0.115672 0.116358 0.11704 0.11772 0.118398 0.119075 0.119752 0.120429 0.121108 0.121789 0.122473 0.12316 0.123851 0.124548 0.125249 0.125957 0.126672 0.127394 0.128123 0.128861 0.129608 0.130364 0.13113 0.131907 0.132694 0.133492 0.134302 0.135124 0.135958 0.136804 0.137664 0.138537 0.139423 0.140323 0.141237 0.142165 0.143108 0.144065 0.145038 0.146024 0.147026 0.148043 0.149076 0.150123 0.151186 0.152264 0.153357 0.154465 0.155589 0.156728 0.157882 0.159051 0.160235 0.161433 0.162646 0.163874 0.165115 0.166371 0.16764 0.168922 0.170218 0.171527 0.172848 0.174181 0.175527 0.176883 0.178251 0.179629 0.181017 0.182415 0.183823 0.185239 0.186663 0.188094 0.189533 0.190978 0.192429 0.193885 0.195346 0.196811 0.198278 0.199748 0.20122 0.202693 0.204166 0.205638 0.207109 0.208577 0.210042 0.211503 0.212959 0.214408 0.215851 0.217286 0.218712 0.220128 0.221533 0.222925 0.224305 0.22567 0.227019 0.228352 0.229667 0.230963 0.232238 0.233492 0.234723 0.235929 0.23711 0.238264 0.23939 0.240486 0.241551 0.242583 0.243581 0.244543 0.245468 0.246354 0.2472 0.248004 0.248765 0.24948 0.250148 0.250768 0.251337 0.251854 0.252317 0.252725 0.253074 0.253365 0.253594 0.25376 0.25386 0.253894];
parallel_jacobi_2_x = [0 0.00440778 0.00871651 0.0129272 0.0170407 0.021058 0.0249801 0.028808 0.0325427 0.0361851 0.0397363 0.0431972 0.0465689 0.0498524 0.0530489 0.0561593 0.0591847 0.0621263 0.0649851 0.0677624 0.0704592 0.0730768 0.0756164 0.0780791 0.0804663 0.0827791 0.0850189 0.0871869 0.0892845 0.091313 0.0932738 0.0951681 0.0969976 0.0987634 0.100467 0.10211 0.103694 0.10522 0.106689 0.108104 0.109465 0.110774 0.112033 0.113243 0.114405 0.115522 0.116594 0.117623 0.118611 0.119558 0.120467 0.121339 0.122175 0.122977 0.123746 0.124483 0.125191 0.125869 0.12652 0.127144 0.127743 0.128319 0.128872 0.129403 0.129913 0.130405 0.130877 0.131333 0.131771 0.132195 0.132603 0.132997 0.133378 0.133747 0.134104 0.134449 0.134784 0.135108 0.135423 0.135728 0.136025 0.136312 0.136591 0.136862 0.137125 0.137381 0.137628 0.137867 0.138099 0.138323 0.138539 0.138746 0.138946 0.139137 0.139319 0.139492 0.139656 0.139811 0.139955 0.140089 0.140212 0.140323 0.140423 0.14051 0.140585 0.140647 0.140694 0.140727 0.140745 0.140748 0.140734 0.140704 0.140656 0.14059 0.140506 0.140402 0.140279 0.140136 0.139971 0.139785 0.139577 0.139346 0.139092 0.138814 0.138512 0.138184 0.137832 0.137453 0.137048 0.136615 0.136155 0.135668 0.135151 0.134606 0.134031 0.133426 0.132791 0.132126 0.131428 0.1307 0.129939 0.129145 0.128319 0.127459 0.126566 0.125638 0.124675 0.123677 0.122643 0.121574 0.120467 0.119324 0.118142 0.116923 0.115665 0.114368 0.11303 0.111652 0.110233 0.108773 0.107269 0.105723 0.104132 0.102497 0.100816 0.0990881 0.0973133 0.0954902 0.0936178 0.0916952 0.0897213 0.0876951 0.0856154 0.0834811 0.081291 0.0790439 0.0767384 0.0743733 0.0719473 0.0694589 0.0669068 0.0642894 0.0616054 0.058853 0.0560309 0.0531374 0.0501708 0.0471296 0.044012 0.0408163 0.0375408 0.0341837 0.0307432 0.0272177 0.0236051 0.0199038 0.0161117 0.0122272 0.00824836 0.00417324 0];
parallel_jacobi_4_x = [0 0.00440778 0.00871651 0.0129272 0.0170407 0.021058 0.0249801 0.028808 0.0325427 0.0361851 0.0397363 0.0431972 0.0465689 0.0498524 0.0530489 0.0561593 0.0591847 0.0621263 0.0649851 0.0677624 0.0704592 0.0730768 0.0756164 0.0780791 0.0804663 0.0827791 0.0850189 0.0871869 0.0892845 0.091313 0.0932738 0.0951681 0.0969976 0.0987634 0.100467 0.10211 0.103694 0.10522 0.106689 0.108104 0.109465 0.110774 0.112033 0.113243 0.114405 0.115522 0.116594 0.117623 0.118611 0.119558 0.120467 0.121339 0.122175 0.122977 0.123746 0.124483 0.125191 0.125869 0.12652 0.127144 0.127743 0.128319 0.128872 0.129403 0.129913 0.130405 0.130877 0.131333 0.131771 0.132195 0.132603 0.132997 0.133378 0.133747 0.134104 0.134449 0.134784 0.135108 0.135423 0.135728 0.136025 0.136312 0.136591 0.136862 0.137125 0.137381 0.137628 0.137867 0.138099 0.138323 0.138539 0.138746 0.138946 0.139137 0.139319 0.139492 0.139656 0.139811 0.139955 0.140089 0.140212 0.140323 0.140423 0.14051 0.140585 0.140647 0.140694 0.140727 0.140745 0.140748 0.140734 0.140704 0.140656 0.14059 0.140506 0.140402 0.140279 0.140136 0.139971 0.139785 0.139577 0.139346 0.139092 0.138814 0.138512 0.138184 0.137832 0.137453 0.137048 0.136615 0.136155 0.135668 0.135151 0.134606 0.134031 0.133426 0.132791 0.132126 0.131428 0.1307 0.129939 0.129145 0.128319 0.127459 0.126566 0.125638 0.124675 0.123677 0.122643 0.121574 0.120467 0.119324 0.118142 0.116923 0.115665 0.114368 0.11303 0.111652 0.110233 0.108773 0.107269 0.105723 0.104132 0.102497 0.100816 0.0990881 0.0973133 0.0954902 0.0936178 0.0916952 0.0897213 0.0876951 0.0856154 0.0834811 0.081291 0.0790439 0.0767384 0.0743733 0.0719473 0.0694589 0.0669068 0.0642894 0.0616054 0.058853 0.0560309 0.0531374 0.0501708 0.0471296 0.044012 0.0408163 0.0375408 0.0341837 0.0307432 0.0272177 0.0236051 0.0199038 0.0161117 0.0122272 0.00824836 0.00417324 0];
parallel_jacobi_4_y = [0 0.00380064 0.00750302 0.0111089 0.0146199 0.0180379 0.0213646 0.0246017 0.0277509 0.030814 0.0337927 0.0366887 0.0395038 0.0422396 0.0448978 0.0474801 0.0499883 0.0524239 0.0547888 0.0570845 0.0593127 0.0614751 0.0635733 0.0656089 0.0675836 0.0694989 0.0713565 0.0731579 0.0749048 0.0765987 0.0782411 0.0798337 0.0813778 0.0828751 0.084327 0.085735 0.0871006 0.0884252 0.0897104 0.0909574 0.0921678 0.0933429 0.0944841 0.0955927 0.0966702 0.0977178 0.0987368 0.0997286 0.100694 0.101635 0.102553 0.103448 0.104323 0.105177 0.106013 0.106831 0.107633 0.108419 0.109191 0.10995 0.110697 0.111433 0.112159 0.112876 0.113585 0.114286 0.114982 0.115672 0.116358 0.11704 0.11772 0.118398 0.119075 0.119752 0.120429 0.121108 0.121789 0.122473 0.12316 0.123851 0.124548 0.125249 0.125957 0.126672 0.127394 0.128123 0.128861 0.129608 0.130364 0.13113 0.131907 0.132694 0.133492 0.134302 0.135124 0.135958 0.136804 0.137664 0.138537 0.139423 0.140323 0.141237 0.142165 0.143108 0.144065 0.145038 0.146024 0.147026 0.148043 0.149076 0.150123 0.151186 0.152264 0.153357 0.154465 0.155589 0.156728 0.157882 0.159051 0.160235 0.161433 0.162646 0.163874 0.165115 0.166371 0.16764 0.168922 0.170218 0.171527 0.172848 0.174181 0.175527 0.176883 0.178251 0.179629 0.181017 0.182415 0.183823 0.185239 0.186663 0.188094 0.189533 0.190978 0.192429 0.193885 0.195346 0.196811 0.198278 0.199748 0.20122 0.202693 0.204166 0.205638 0.207109 0.208577 0.210042 0.211503 0.212959 0.214408 0.215851 0.217286 0.218712 0.220128 0.221533 0.222925 0.224305 0.22567 0.227019 0.228352 0.229667 0.230963 0.232238 0.233492 0.234723 0.235929 0.23711 0.238264 0.23939 0.240486 0.241551 0.242583 0.243581 0.244543 0.245468 0.246354 0.2472 0.248004 0.248765 0.24948 0.250148 0.250768 0.251337 0.251854 0.252317 0.252725 0.253074 0.253365 0.253594 0.25376 0.25386 0.253894];
parallel_jacobi_8_x = [0 0.00440778 0.00871651 0.0129272 0.0170407 0.021058 0.0249801 0.028808 0.0325427 0.0361851 0.0397363 0.0431972 0.0465689 0.0498524 0.0530489 0.0561593 0.0591847 0.0621263 0.0649851 0.0677624 0.0704592 0.0730768 0.0756164 0.0780791 0.0804663 0.0827791 0.0850189 0.0871869 0.0892845 0.091313 0.0932738 0.0951681 0.0969976 0.0987634 0.100467 0.10211 0.103694 0.10522 0.106689 0.108104 0.109465 0.110774 0.112033 0.113243 0.114405 0.115522 0.116594 0.117623 0.118611 0.119558 0.120467 0.121339 0.122175 0.122977 0.123746 0.124483 0.125191 0.125869 0.12652 0.127144 0.127743 0.128319 0.128872 0.129403 0.129913 0.130405 0.130877 0.131333 0.131771 0.132195 0.132603 0.132997 0.133378 0.133747 0.134104 0.134449 0.134784 0.135108 0.135423 0.135728 0.136025 0.136312 0.136591 0.136862 0.137125 0.137381 0.137628 0.137867 0.138099 0.138323 0.138539 0.138746 0.138946 0.139137 0.139319 0.139492 0.139656 0.139811 0.139955 0.140089 0.140212 0.140323 0.140423 0.14051 0.140585 0.140647 0.140694 0.140727 0.140745 0.140748 0.140734 0.140704 0.140656 0.14059 0.140506 0.140402 0.140279 0.140136 0.139971 0.139785 0.139577 0.139346 0.139092 0.138814 0.138512 0.138184 0.137832 0.137453 0.137048 0.136615 0.136155 0.135668 0.135151 0.134606 0.134031 0.133426 0.132791 0.132126 0.131428 0.1307 0.129939 0.129145 0.128319 0.127459 0.126566 0.125638 0.124675 0.123677 0.122643 0.121574 0.120467 0.119324 0.118142 0.116923 0.115665 0.114368 0.11303 0.111652 0.110233 0.108773 0.107269 0.105723 0.104132 0.102497 0.100816 0.0990881 0.0973133 0.0954902 0.0936178 0.0916952 0.0897213 0.0876951 0.0856154 0.0834811 0.081291 0.0790439 0.0767384 0.0743733 0.0719473 0.0694589 0.0669068 0.0642894 0.0616054 0.058853 0.0560309 0.0531374 0.0501708 0.0471296 0.044012 0.0408163 0.0375408 0.0341837 0.0307432 0.0272177 0.0236051 0.0199038 0.0161117 0.0122272 0.00824836 0.00417324 0];
parallel_jacobi_8_y = [0 0.00380064 0.00750302 0.0111089 0.0146199 0.0180379 0.0213646 0.0246017 0.0277509 0.030814 0.0337927 0.0366887 0.0395038 0.0422396 0.0448978 0.0474801 0.0499883 0.0524239 0.0547888 0.0570845 0.0593127 0.0614751 0.0635733 0.0656089 0.0675836 0.0694989 0.0713565 0.0731579 0.0749048 0.0765987 0.0782411 0.0798337 0.0813778 0.0828751 0.084327 0.085735 0.0871006 0.0884252 0.0897104 0.0909574 0.0921678 0.0933429 0.0944841 0.0955927 0.0966702 0.0977178 0.0987368 0.0997286 0.100694 0.101635 0.102553 0.103448 0.104323 0.105177 0.106013 0.106831 0.107633 0.108419 0.109191 0.10995 0.110697 0.111433 0.112159 0.112876 0.113585 0.114286 0.114982 0.115672 0.116358 0.11704 0.11772 0.118398 0.119075 0.119752 0.120429 0.121108 0.121789 0.122473 0.12316 0.123851 0.124548 0.125249 0.125957 0.126672 0.127394 0.128123 0.128861 0.129608 0.130364 0.13113 0.131907 0.132694 0.133492 0.134302 0.135124 0.135958 0.136804 0.137664 0.138537 0.139423 0.140323 0.141237 0.142165 0.143108 0.144065 0.145038 0.146024 0.147026 0.148043 0.149076 0.150123 0.151186 0.152264 0.153357 0.154465 0.155589 0.156728 0.157882 0.159051 0.160235 0.161433 0.162646 0.163874 0.165115 0.166371 0.16764 0.168922 0.170218 0.171527 0.172848 0.174181 0.175527 0.176883 0.178251 0.179629 0.181017 0.182415 0.183823 0.185239 0.186663 0.188094 0.189533 0.190978 0.192429 0.193885 0.195346 0.196811 0.198278 0.199748 0.20122 0.202693 0.204166 0.205638 0.207109 0.208577 0.210042 0.211503 0.212959 0.214408 0.215851 0.217286 0.218712 0.220128 0.221533 0.222925 0.224305 0.22567 0.227019 0.228352 0.229667 0.230963 0.232238 0.233492 0.234723 0.235929 0.23711 0.238264 0.23939 0.240486 0.241551 0.242583 0.243581 0.244543 0.245468 0.246354 0.2472 0.248004 0.248765 0.24948 0.250148 0.250768 0.251337 0.251854 0.252317 0.252725 0.253074 0.253365 0.253594 0.25376 0.25386 0.253894];
x = linspace(-1,1,201);

figure
subplot(1,2,1)
plot(x, parallel_jacobi_2_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_jacobi_1_y, 'b--', 'LineWidth', 1);
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Jacobi Parallel Solution with 2 processors at y = 0');
legend('Parallel Solution', 'Serial Solution');

subplot(1,2,2)
plot(x, parallel_jacobi_2_x , 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_jacobi_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Jacobi Parallel Solution with 2 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');

figure
subplot(1,2,1)
plot(x, parallel_jacobi_4_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_jacobi_1_y, 'b--', 'LineWidth', 1);
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Jacobi Parallel Solution with 4 processors at y = 0');
legend('Parallel Solution', 'Serial Solution');

subplot(1,2,2)
plot(x, parallel_jacobi_4_x, 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_jacobi_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Jacobi Parallel Solution with 4 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');

figure
subplot(1,2,1)
plot(x, parallel_jacobi_8_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_jacobi_1_y, 'b--', 'LineWidth', 1); 
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Jacobi Parallel Solution with 8 processors at y = 0');
legend('Parallel Solution', 'Serial Solution')

subplot(1,2,2)
plot(x, parallel_jacobi_8_x, 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_jacobi_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Jacobi Parallel Solution with 8 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');