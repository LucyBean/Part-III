import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage

def emailText(subject, body):
    msg = MIMEText(body)
        
    me = you = "lucymair@live.co.uk"
    
    msg["Subject"] = subject
    msg["From"] = me
    msg["To"] = you
    
    server = smtplib.SMTP("smtp-mail.outlook.com", 587)
    server.ehlo()
    server.starttls()
    server.login("lucymair@live.co.uk","uukulxcxvvisnfrw")
    server.sendmail(me, [you], msg.as_string())

def emailImage(subject, body, imagePath):
    msg = MIMEMultipart()
    me = you = "lucymair@live.co.uk"
    
    msg["Subject"] = subject
    msg["From"] = me
    msg["To"] = you
    
    text = MIMEText(body)
    msg.attach(text)
    
    with open(imagePath, "rb") as f:
        img = MIMEImage(f.read())
        msg.attach(img)

    server = smtplib.SMTP("smtp-mail.outlook.com", 587)
    server.ehlo()
    server.starttls()
    server.login("lucymair@live.co.uk","uukulxcxvvisnfrw")
    server.sendmail(me, [you], msg.as_string())
