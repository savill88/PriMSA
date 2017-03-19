from django.conf.urls import url
from primsa import views

urlpatterns=[
    url(r'^$',views.index, name='index'),
    url(r'^ncbi/$',views.ncbi, name='ncbi'),
    url(r'^msa/$',views.msa, name='msa'),
]
