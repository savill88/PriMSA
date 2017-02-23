from django.conf.urls import url
from primsa import views

urlpatterns=[
    url(r'^$',views.index, name='index'),
]
