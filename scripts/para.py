#!/usr/bin/env python3
"""para.py - paracloud SSH helper (password from ~/.ssh/paracloud_pass).
Usage: python para.py run "<cmd>"   |   python para.py put <local> <remote>"""
import sys, paramiko
HOST, PORT, USER = "ssh.paracloud.com", 2222, "pxyl058@GUANGZHOUXYL"
pw = open(r"C:\Users\13614\.ssh\paracloud_pass").read()
c = paramiko.SSHClient(); c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
c.connect(HOST, port=PORT, username=USER, password=pw, timeout=25, allow_agent=False, look_for_keys=False)
if sys.argv[1] == "run":
    _, out, err = c.exec_command(sys.argv[2], timeout=1800)
    sys.stdout.write(out.read().decode(errors="replace"))
    e = err.read().decode(errors="replace")
    if e.strip(): sys.stderr.write("[stderr] " + e[-2000:])
    sys.exit(out.channel.recv_exit_status())
elif sys.argv[1] == "put":
    s = c.open_sftp(); s.put(sys.argv[2], sys.argv[3]); s.close(); print("uploaded", sys.argv[3])
c.close()
