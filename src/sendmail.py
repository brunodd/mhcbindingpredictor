import smtplib

def alert(content):
    print("Sending mail:", content)
    fromaddr = 'dedekenbruno@gmail.com'
    toaddrs  = 'dedekenbruno@gmail.com'
    subject = 'Training alert'
    username = fromaddr
    password = 'MKt^R2q1'
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.ehlo()
    server.starttls()
    server.login(username,password)
    msg = "\r\n".join([
      "From: " + fromaddr,
      "To: " + toaddrs,
      "Subject: " + subject,
      "",
      content
      ])
    server.sendmail(fromaddr, toaddrs, msg)
    server.quit()
    print("Mail sent")

if __name__ == '__main__':
    alert("Some alert")
