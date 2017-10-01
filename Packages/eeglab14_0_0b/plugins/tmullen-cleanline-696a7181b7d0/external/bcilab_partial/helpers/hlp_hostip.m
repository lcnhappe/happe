function hostip = hlp_hostip()
% Find the host IP address in a platform-independent fashion.
import java.net.InetAddress;
hostip = char(InetAddress.getLocalHost());
hostip = hostip(find(hostip=='/',1)+1:end);