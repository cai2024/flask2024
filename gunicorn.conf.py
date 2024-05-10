workers = 8
bind = '0.0.0.0:6090'
accesslog = './app_logs/gunicorn_access.log'
errorlog = './app_logs/gunicorn_error.log'
loglevel = 'warning'
backlog = 2048
