function hostname = hlp_hostname()
% Find the hostname in a platform-independent fashion.
import java.net.InetAddress;
hostname = char(InetAddress.getLocalHost());
hostname = hostname(1:find(hostname=='/',1)-1);