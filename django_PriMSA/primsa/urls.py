from django.conf.urls import url
from primsa import views

urlpatterns=[
    url(r'^$',views.index, name='index'),
    url(r'^retrieve/$',views.retrieve, name='retrieve'),
    url(r'^msa/$',views.msa, name='msa'),
]
