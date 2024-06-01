clc
clear
serial_gs_1_x = [0 0.00443216 0.00876529 0.0130004 0.0171384 0.0211803 0.025127 0.0289796 0.0327391 0.0364064 0.0399827 0.0434688 0.0468658 0.0501749 0.053397 0.0565333 0.0595848 0.0625526 0.065438 0.0682421 0.0709659 0.0736108 0.0761779 0.0786684 0.0810836 0.0834248 0.0856933 0.0878903 0.0900172 0.0920754 0.0940661 0.0959907 0.0978508 0.0996475 0.101382 0.103057 0.104673 0.106231 0.107733 0.10918 0.110575 0.111918 0.11321 0.114454 0.115651 0.116802 0.117909 0.118974 0.119997 0.12098 0.121924 0.122832 0.123704 0.124542 0.125347 0.126121 0.126864 0.127579 0.128266 0.128926 0.129562 0.130173 0.130762 0.131328 0.131874 0.1324 0.132907 0.133397 0.133869 0.134326 0.134767 0.135194 0.135607 0.136006 0.136393 0.136768 0.137132 0.137485 0.137827 0.138158 0.13848 0.138792 0.139095 0.139389 0.139673 0.139949 0.140215 0.140473 0.140721 0.140961 0.141191 0.141412 0.141624 0.141825 0.142017 0.142198 0.142369 0.142529 0.142677 0.142813 0.142937 0.143049 0.143147 0.143231 0.143302 0.143358 0.143398 0.143423 0.143431 0.143423 0.143397 0.143353 0.14329 0.143208 0.143107 0.142985 0.142842 0.142677 0.142491 0.142282 0.14205 0.141794 0.141514 0.14121 0.14088 0.140524 0.140142 0.139733 0.139297 0.138833 0.138341 0.137821 0.137272 0.136693 0.136084 0.135445 0.134775 0.134074 0.133341 0.132577 0.13178 0.13095 0.130087 0.129191 0.128261 0.127297 0.126297 0.125263 0.124193 0.123087 0.121944 0.120764 0.119547 0.118292 0.116999 0.115666 0.114293 0.112881 0.111427 0.109932 0.108395 0.106814 0.10519 0.103522 0.101808 0.100048 0.0982412 0.0963863 0.0944826 0.0925289 0.0905243 0.0884677 0.086358 0.0841939 0.0819744 0.0796982 0.077364 0.0749704 0.0725162 0.07 0.0674202 0.0647755 0.0620643 0.0592852 0.0564364 0.0535165 0.0505237 0.0474564 0.044313 0.0410916 0.0377906 0.0344082 0.0309425 0.0273918 0.0237542 0.0200279 0.016211 0.0123016 0.00829793 0.00419802 0];
serial_gs_1_y = [0 0.00380793 0.00751763 0.0111309 0.0146494 0.018075 0.0214095 0.0246546 0.027812 0.0308835 0.0338708 0.0367758 0.0396001 0.0423454 0.0450136 0.0476062 0.0501251 0.0525719 0.0549483 0.057256 0.0594967 0.0616721 0.0637838 0.0658334 0.0678225 0.0697529 0.0716262 0.0734438 0.0752074 0.0769186 0.078579 0.08019 0.0817532 0.0832702 0.0847424 0.0861713 0.0875584 0.0889051 0.0902129 0.0914832 0.0927175 0.093917 0.0950832 0.0962175 0.0973211 0.0983954 0.0994417 0.100461 0.101455 0.102425 0.103372 0.104297 0.105202 0.106087 0.106954 0.107804 0.108638 0.109457 0.110262 0.111054 0.111835 0.112605 0.113365 0.114117 0.11486 0.115597 0.116328 0.117054 0.117776 0.118495 0.119211 0.119926 0.12064 0.121354 0.122068 0.122785 0.123503 0.124224 0.124949 0.125678 0.126412 0.127152 0.127898 0.12865 0.12941 0.130177 0.130953 0.131738 0.132532 0.133336 0.13415 0.134975 0.13581 0.136658 0.137517 0.138388 0.139272 0.140168 0.141078 0.142001 0.142937 0.143888 0.144852 0.145831 0.146823 0.147831 0.148853 0.14989 0.150942 0.152009 0.15309 0.154187 0.155299 0.156426 0.157567 0.158724 0.159896 0.161082 0.162283 0.163498 0.164728 0.165973 0.167231 0.168503 0.169789 0.171088 0.1724 0.173725 0.175063 0.176412 0.177774 0.179147 0.180531 0.181926 0.183331 0.184746 0.18617 0.187603 0.189044 0.190493 0.191949 0.193412 0.194881 0.196356 0.197835 0.199318 0.200805 0.202295 0.203787 0.20528 0.206773 0.208266 0.209758 0.211249 0.212736 0.21422 0.215699 0.217173 0.218641 0.220101 0.221552 0.222995 0.224426 0.225847 0.227255 0.228649 0.230028 0.231392 0.232738 0.234066 0.235375 0.236662 0.237928 0.23917 0.240388 0.241579 0.242744 0.243879 0.244985 0.246058 0.247099 0.248105 0.249075 0.250007 0.250901 0.251753 0.252563 0.253329 0.25405 0.254723 0.255347 0.25592 0.25644 0.256907 0.257317 0.257669 0.257961 0.258192 0.258359 0.25846 0.258494];
parallel_gs_2_y = [0 0.0038405 0.00758273 0.0112284 0.0147793 0.0182371 0.0216036 0.0248805 0.0280696 0.0311724 0.0341909 0.0371266 0.0399813 0.0427567 0.0454546 0.0480765 0.0506242 0.0530994 0.0555038 0.0578388 0.0601065 0.0623081 0.0644456 0.0665204 0.0685343 0.0704887 0.0723854 0.0742257 0.0760116 0.0777442 0.0794255 0.0810566 0.0826395 0.0841752 0.0856656 0.0871119 0.0885159 0.0898786 0.0912019 0.0924869 0.0937353 0.0949481 0.0961271 0.0972732 0.0983883 0.0994731 0.10053 0.101558 0.102561 0.103539 0.104494 0.105426 0.106337 0.107227 0.108099 0.108954 0.109791 0.110614 0.111422 0.112216 0.112999 0.11377 0.114531 0.115283 0.116027 0.116763 0.117493 0.118218 0.118938 0.119654 0.120368 0.12108 0.12179 0.122501 0.123211 0.123923 0.124637 0.125353 0.126073 0.126797 0.127525 0.128259 0.128999 0.129745 0.130498 0.131259 0.132028 0.132805 0.133592 0.134388 0.135195 0.136012 0.13684 0.137679 0.138531 0.139394 0.14027 0.141158 0.142059 0.142973 0.143902 0.144843 0.145799 0.146769 0.147754 0.148752 0.149766 0.150794 0.151838 0.152896 0.153969 0.155057 0.15616 0.157278 0.158412 0.159559 0.160723 0.1619 0.163093 0.1643 0.165522 0.166758 0.168008 0.169272 0.17055 0.171841 0.173145 0.174462 0.175792 0.177134 0.178488 0.179853 0.18123 0.182617 0.184015 0.185422 0.186839 0.188265 0.189699 0.191141 0.192591 0.194047 0.195509 0.196977 0.19845 0.199927 0.201408 0.202891 0.204377 0.205863 0.207351 0.208838 0.210325 0.211809 0.213292 0.21477 0.216244 0.217712 0.219175 0.220629 0.222077 0.223514 0.224941 0.226357 0.22776 0.22915 0.230525 0.231884 0.233227 0.23455 0.235855 0.237139 0.238401 0.239639 0.240854 0.242042 0.243203 0.244335 0.245438 0.246508 0.247546 0.248549 0.249517 0.250446 0.251337 0.252187 0.252995 0.253759 0.254477 0.255148 0.255771 0.256342 0.256861 0.257325 0.257734 0.258085 0.258376 0.258605 0.258771 0.258871 0.258905];
parallel_gs_2_x = [0 0.004464 0.00882892 0.0130958 0.0172654 0.0213388 0.0253169 0.0292008 0.0329912 0.0366894 0.0402961 0.0438125 0.0472395 0.0505782 0.0538296 0.0569948 0.0600747 0.0630707 0.0659836 0.0688148 0.0715651 0.0742361 0.0768287 0.0793442 0.0817837 0.0841487 0.0864402 0.0886597 0.0908083 0.0928876 0.0948985 0.0968429 0.0987216 0.100537 0.102289 0.10398 0.105611 0.107184 0.108701 0.110162 0.111568 0.112923 0.114226 0.115481 0.116687 0.117846 0.118961 0.120032 0.12106 0.122049 0.122998 0.123909 0.124784 0.125624 0.12643 0.127205 0.127948 0.128662 0.129347 0.130006 0.130639 0.131247 0.131831 0.132394 0.132935 0.133456 0.133958 0.134441 0.134908 0.135358 0.135792 0.136211 0.136617 0.137009 0.137388 0.137756 0.138111 0.138456 0.13879 0.139114 0.139428 0.139733 0.140028 0.140315 0.140592 0.140861 0.14112 0.141372 0.141614 0.141848 0.142073 0.14229 0.142496 0.142694 0.142882 0.14306 0.143228 0.143386 0.143532 0.143667 0.14379 0.143902 0.144 0.144086 0.144157 0.144215 0.144257 0.144285 0.144296 0.144291 0.144269 0.144229 0.144171 0.144095 0.143998 0.143883 0.143746 0.143588 0.143408 0.143206 0.142981 0.142732 0.142459 0.142162 0.141839 0.141491 0.141116 0.140715 0.140286 0.139829 0.139344 0.138831 0.138287 0.137715 0.137111 0.136478 0.135813 0.135117 0.134388 0.133628 0.132834 0.132008 0.131147 0.130252 0.129323 0.128359 0.127359 0.126324 0.125252 0.124144 0.122998 0.121814 0.120592 0.119332 0.118032 0.116692 0.115311 0.11389 0.112427 0.110922 0.109373 0.107781 0.106144 0.104462 0.102733 0.100958 0.099135 0.0972634 0.0953419 0.09337 0.0913461 0.0892697 0.0871392 0.0849539 0.0827123 0.0804135 0.0780558 0.0756384 0.0731596 0.0706183 0.0680127 0.0653418 0.0626038 0.0597975 0.0569209 0.0539729 0.0509515 0.0478554 0.0446825 0.0414316 0.0381006 0.034688 0.0311918 0.0276104 0.0239419 0.0201846 0.0163365 0.0123958 0.00836078 0.00422946 0];
parallel_gs_4_x = [0 0.004464 0.00882892 0.0130958 0.0172654 0.0213388 0.0253169 0.0292008 0.0329912 0.0366894 0.0402961 0.0438125 0.0472395 0.0505782 0.0538296 0.0569948 0.0600747 0.0630707 0.0659836 0.0688148 0.0715651 0.0742361 0.0768287 0.0793442 0.0817837 0.0841487 0.0864402 0.0886597 0.0908083 0.0928876 0.0948985 0.0968429 0.0987216 0.100537 0.102289 0.10398 0.105611 0.107184 0.108701 0.110162 0.111568 0.112923 0.114226 0.115481 0.116687 0.117846 0.118961 0.120032 0.12106 0.122049 0.122998 0.123909 0.124784 0.125624 0.12643 0.127205 0.127948 0.128662 0.129347 0.130006 0.130639 0.131247 0.131831 0.132394 0.132935 0.133456 0.133958 0.134441 0.134908 0.135358 0.135792 0.136211 0.136617 0.137009 0.137388 0.137756 0.138111 0.138456 0.13879 0.139114 0.139428 0.139733 0.140028 0.140315 0.140592 0.140861 0.14112 0.141372 0.141614 0.141848 0.142073 0.14229 0.142496 0.142694 0.142882 0.14306 0.143228 0.143386 0.143532 0.143667 0.14379 0.143902 0.144 0.144086 0.144157 0.144215 0.144257 0.144285 0.144296 0.144291 0.144269 0.144229 0.144171 0.144095 0.143998 0.143883 0.143746 0.143588 0.143408 0.143206 0.142981 0.142732 0.142459 0.142162 0.141839 0.141491 0.141116 0.140715 0.140286 0.139829 0.139344 0.138831 0.138287 0.137715 0.137111 0.136478 0.135813 0.135117 0.134388 0.133628 0.132834 0.132008 0.131147 0.130252 0.129323 0.128359 0.127359 0.126324 0.125252 0.124144 0.122998 0.121814 0.120592 0.119332 0.118032 0.116692 0.115311 0.11389 0.112427 0.110922 0.109373 0.107781 0.106144 0.104462 0.102733 0.100958 0.099135 0.0972634 0.0953419 0.09337 0.0913461 0.0892697 0.0871392 0.0849539 0.0827123 0.0804135 0.0780558 0.0756384 0.0731596 0.0706183 0.0680127 0.0653418 0.0626038 0.0597975 0.0569209 0.0539729 0.0509515 0.0478554 0.0446825 0.0414316 0.0381006 0.034688 0.0311918 0.0276104 0.0239419 0.0201846 0.0163365 0.0123958 0.00836078 0.00422946 0];
parallel_gs_4_y = [0 0.0038405 0.00758273 0.0112284 0.0147793 0.0182371 0.0216036 0.0248805 0.0280696 0.0311724 0.0341909 0.0371266 0.0399813 0.0427567 0.0454546 0.0480765 0.0506242 0.0530994 0.0555038 0.0578388 0.0601065 0.0623081 0.0644456 0.0665204 0.0685343 0.0704887 0.0723854 0.0742257 0.0760116 0.0777442 0.0794255 0.0810566 0.0826395 0.0841752 0.0856656 0.0871119 0.0885159 0.0898786 0.0912019 0.0924869 0.0937353 0.0949481 0.0961271 0.0972732 0.0983883 0.0994731 0.10053 0.101558 0.102561 0.103539 0.104494 0.105426 0.106337 0.107227 0.108099 0.108954 0.109791 0.110614 0.111422 0.112216 0.112999 0.11377 0.114531 0.115283 0.116027 0.116763 0.117493 0.118218 0.118938 0.119654 0.120368 0.12108 0.12179 0.122501 0.123211 0.123923 0.124637 0.125353 0.126073 0.126797 0.127525 0.128259 0.128999 0.129745 0.130498 0.131259 0.132028 0.132805 0.133592 0.134388 0.135195 0.136012 0.13684 0.137679 0.138531 0.139394 0.14027 0.141158 0.142059 0.142973 0.143902 0.144843 0.145799 0.146769 0.147754 0.148752 0.149766 0.150794 0.151838 0.152896 0.153969 0.155057 0.15616 0.157278 0.158412 0.159559 0.160723 0.1619 0.163093 0.1643 0.165522 0.166758 0.168008 0.169272 0.17055 0.171841 0.173145 0.174462 0.175792 0.177134 0.178488 0.179853 0.18123 0.182617 0.184015 0.185422 0.186839 0.188265 0.189699 0.191141 0.192591 0.194047 0.195509 0.196977 0.19845 0.199927 0.201408 0.202891 0.204377 0.205863 0.207351 0.208838 0.210325 0.211809 0.213292 0.21477 0.216244 0.217712 0.219175 0.220629 0.222077 0.223514 0.224941 0.226357 0.22776 0.22915 0.230525 0.231884 0.233227 0.23455 0.235855 0.237139 0.238401 0.239639 0.240854 0.242042 0.243203 0.244335 0.245438 0.246508 0.247546 0.248549 0.249517 0.250446 0.251337 0.252187 0.252995 0.253759 0.254477 0.255148 0.255771 0.256342 0.256861 0.257325 0.257734 0.258085 0.258376 0.258605 0.258771 0.258871 0.258905];
parallel_gs_8_x = [0 0.004464 0.00882892 0.0130958 0.0172654 0.0213388 0.0253169 0.0292008 0.0329912 0.0366894 0.0402961 0.0438125 0.0472395 0.0505782 0.0538296 0.0569948 0.0600747 0.0630707 0.0659836 0.0688148 0.0715651 0.0742361 0.0768287 0.0793442 0.0817837 0.0841487 0.0864402 0.0886597 0.0908083 0.0928876 0.0948985 0.0968429 0.0987216 0.100537 0.102289 0.10398 0.105611 0.107184 0.108701 0.110162 0.111568 0.112923 0.114226 0.115481 0.116687 0.117846 0.118961 0.120032 0.12106 0.122049 0.122998 0.123909 0.124784 0.125624 0.12643 0.127205 0.127948 0.128662 0.129347 0.130006 0.130639 0.131247 0.131831 0.132394 0.132935 0.133456 0.133958 0.134441 0.134908 0.135358 0.135792 0.136211 0.136617 0.137009 0.137388 0.137756 0.138111 0.138456 0.13879 0.139114 0.139428 0.139733 0.140028 0.140315 0.140592 0.140861 0.14112 0.141372 0.141614 0.141848 0.142073 0.14229 0.142496 0.142694 0.142882 0.14306 0.143228 0.143386 0.143532 0.143667 0.14379 0.143902 0.144 0.144086 0.144157 0.144215 0.144257 0.144285 0.144296 0.144291 0.144269 0.144229 0.144171 0.144095 0.143998 0.143883 0.143746 0.143588 0.143408 0.143206 0.142981 0.142732 0.142459 0.142162 0.141839 0.141491 0.141116 0.140715 0.140286 0.139829 0.139344 0.138831 0.138287 0.137715 0.137111 0.136478 0.135813 0.135117 0.134388 0.133628 0.132834 0.132008 0.131147 0.130252 0.129323 0.128359 0.127359 0.126324 0.125252 0.124144 0.122998 0.121814 0.120592 0.119332 0.118032 0.116692 0.115311 0.11389 0.112427 0.110922 0.109373 0.107781 0.106144 0.104462 0.102733 0.100958 0.099135 0.0972634 0.0953419 0.09337 0.0913461 0.0892697 0.0871392 0.0849539 0.0827123 0.0804135 0.0780558 0.0756384 0.0731596 0.0706183 0.0680127 0.0653418 0.0626038 0.0597975 0.0569209 0.0539729 0.0509515 0.0478554 0.0446825 0.0414316 0.0381006 0.034688 0.0311918 0.0276104 0.0239419 0.0201846 0.0163365 0.0123958 0.00836078 0.00422946 0];
parallel_gs_8_y = [0 0.0038405 0.00758273 0.0112284 0.0147793 0.0182371 0.0216036 0.0248805 0.0280696 0.0311724 0.0341909 0.0371266 0.0399813 0.0427567 0.0454546 0.0480765 0.0506242 0.0530994 0.0555038 0.0578388 0.0601065 0.0623081 0.0644456 0.0665204 0.0685343 0.0704887 0.0723854 0.0742257 0.0760116 0.0777442 0.0794255 0.0810566 0.0826395 0.0841752 0.0856656 0.0871119 0.0885159 0.0898786 0.0912019 0.0924869 0.0937353 0.0949481 0.0961271 0.0972732 0.0983883 0.0994731 0.10053 0.101558 0.102561 0.103539 0.104494 0.105426 0.106337 0.107227 0.108099 0.108954 0.109791 0.110614 0.111422 0.112216 0.112999 0.11377 0.114531 0.115283 0.116027 0.116763 0.117493 0.118218 0.118938 0.119654 0.120368 0.12108 0.12179 0.122501 0.123211 0.123923 0.124637 0.125353 0.126073 0.126797 0.127525 0.128259 0.128999 0.129745 0.130498 0.131259 0.132028 0.132805 0.133592 0.134388 0.135195 0.136012 0.13684 0.137679 0.138531 0.139394 0.14027 0.141158 0.142059 0.142973 0.143902 0.144843 0.145799 0.146769 0.147754 0.148752 0.149766 0.150794 0.151838 0.152896 0.153969 0.155057 0.15616 0.157278 0.158412 0.159559 0.160723 0.1619 0.163093 0.1643 0.165522 0.166758 0.168008 0.169272 0.17055 0.171841 0.173145 0.174462 0.175792 0.177134 0.178488 0.179853 0.18123 0.182617 0.184015 0.185422 0.186839 0.188265 0.189699 0.191141 0.192591 0.194047 0.195509 0.196977 0.19845 0.199927 0.201408 0.202891 0.204377 0.205863 0.207351 0.208838 0.210325 0.211809 0.213292 0.21477 0.216244 0.217712 0.219175 0.220629 0.222077 0.223514 0.224941 0.226357 0.22776 0.22915 0.230525 0.231884 0.233227 0.23455 0.235855 0.237139 0.238401 0.239639 0.240854 0.242042 0.243203 0.244335 0.245438 0.246508 0.247546 0.248549 0.249517 0.250446 0.251337 0.252187 0.252995 0.253759 0.254477 0.255148 0.255771 0.256342 0.256861 0.257325 0.257734 0.258085 0.258376 0.258605 0.258771 0.258871 0.258905];
x = linspace(-1,1,201);

figure
subplot(1,2,1)
plot(x, parallel_gs_2_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_gs_1_y, 'b--', 'LineWidth', 1);
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Gauss Seidel Red Black Solution with 2 processors at y = 0');
legend('Parallel Solution', 'Serial Solution');

subplot(1,2,2)
plot(x, parallel_gs_2_x , 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_gs_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Gauss Seidel Red Black Solution with 2 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');

figure
subplot(1,2,1)
plot(x, parallel_gs_4_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_gs_1_y, 'b--', 'LineWidth', 1);
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Gauss Seidel Red Black Solution with 4 processors at y = 0');
legend('Parallel Solution', 'Serial Solution');

subplot(1,2,2)
plot(x, parallel_gs_4_x, 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_gs_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Gauss Seidel Red Black Solution with 4 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');

figure
subplot(1,2,1)
plot(x, parallel_gs_8_y, 'r--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
hold on
plot(x, serial_gs_1_y, 'b--', 'LineWidth', 1); 
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Gauss Seidel Red Black Solution with 8 processors at y = 0');
legend('Parallel Solution', 'Serial Solution')

subplot(1,2,2)
plot(x, parallel_gs_8_x, 'r--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
hold on
plot(x, serial_gs_1_x, 'b--', 'LineWidth', 1);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Gauss Seidel Red Black Solution with 8 processors at x = 0');
legend('Parallel Solution', 'Serial Solution');